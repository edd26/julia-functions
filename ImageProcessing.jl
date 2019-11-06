using Statistics
using Combinatorics
using ImageFiltering

"""
    rotate_img_around_center(img, angle = 5pi/6)

Function rotates a single image (or a frame) around the center of the image by
@angle radians.
"""
function rotate_img_around_center(img, angle = 5pi/6)
  θ = angle
  rot = recenter(RotMatrix(θ), [size(img)...] .÷ 2)  # a rotation around the center
  x_translation = 0
  y_translation = 0
  tform = rot ∘ Translation(y_translation, x_translation)
  img2 = warp(img, rot, axes(img))

  return img2
end


"""
    get_local_gradients(video_array, centers, sub_img_size)

Computes the gradients in the subimage, takes the mean of sum of absolute
values of both hotizontal and vertical gradients as a representative of a
subimage.
"""
function get_local_gradients(video_array, centers, sub_img_size)
    @debug "Entering get_local_gradients"
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size
    h, w, len = get_video_dimension(video_array)
    extracted_pixels = zeros(sub_img_size, sub_img_size, len)

    @debug "starting frame procesing"
    for frame = 1:len
        img = video_array[frame]
        img_grad = imgradients(img, KernelFactors.ando3, "replicate")
        img_grad_abs = map(abs, img_grad[1]) + map(abs, img_grad[2])
        for index_x = 1:size(centers,2)
            c_x = centers[2, index_x]
            for index_y = 1:size(centers,2)
                c_y = centers[1, index_y]
                sub_img  = img_grad_abs[(c_x-half_size):(c_x+half_size),
                                        (c_y-half_size):(c_y+half_size)]

                extracted_pixels[index_x, index_y, frame] =mean(sub_img)
            end
        end
        # @debug "Next frame" frame
    end
    return extracted_pixels
end


"""
    get_img_gradient(img)

Computes the gradients in the `img`.
"""
function get_img_gradient(img)
    @debug "Entering get_local_gradients"

    img_grad = imgradients(img, KernelFactors.ando3, "replicate")

    grad_1 = img_grad[1] .+ abs(findmin(img_grad[1])[1])
    grad_1 ./= findmax(grad_1)[1]

    grad_2 = img_grad[2] .+ abs(findmin(img_grad[2])[1])
    grad_2 ./= findmax(grad_2)[1]

    grad_sum = grad_1 + grad_2
    grad_sum .+= abs(findmin(grad_sum)[1])
    grad_sum ./= findmax(grad_sum)[1]

    return grad_sum
end


"""
    get_local_img_correlations(img, centers, sub_img_size, shift;
                                                        with_gradient=false)

Computes the correlation between the subimages and subimages shifted by values
from range -`shift`:`shift` and returns array of size
length(`centers`) x length(`centers`).

Each of the subimage is center around values stored in  @centers
"""
function get_local_img_correlations(img, centers, sub_img_size;
                                                    with_gradient=false)
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size#
    h, w = size(img)
    extracted_pixels = zeros(sub_img_size, sub_img_size)
    local_correlation = zeros(size(centers,1))

    if with_gradient
        img = get_img_gradient(img)
    end

    position = 1;
    for index = centers
        c_x = index[1]
        c_y = index[2]
        subimage = img[(c_x-half_range):(c_x+half_range),
                        (c_y-half_range):(c_y+half_range)]
        center = img[(c_x-half_size):(c_x+half_size), (c_y-half_size):(c_y+half_size)]

        corelation = center .* subimage
        corelation = sum(corelation)
        local_correlation[position] += corelation

        local_correlation[position] /= 256*(sub_img_size^2)^2
        position += 1;
    end

    return local_correlation
end




"""
    get_local_img_correlations(img,centers, masks)

Takes `img` and computes crosscorrelation with set of `masks` around the
`centers`. Crosscorrelation is computed as convolution of the mask and the area
around coordinates stored in `centres`.
"""
function get_local_img_correlations(img, centers, masks; with_gradient=false)
    masks_num = length(masks)
    sub_img_size = size(masks[1],1)
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size
    h, w = size(img)
    local_correlation = zeros(masks_num, size(centers,1) )
    index = centers[1]
    masks_num = length(masks)

    if with_gradient
        img = get_img_gradient(img)
    end

    position = 1;
    for index = centers
        c_x = index[1]
        c_y = index[2]

        center = img[(c_x-half_size):(c_x+half_size), (c_y-half_size):(c_y+half_size)]
        mask_pos = 1
        for mask in masks

            corelation = center .* mask

            corelation = sum(corelation)
            local_correlation[mask_pos, position] += corelation
            local_correlation[mask_pos, position] /= (sub_img_size^2)
            # local_correlation[position, mask_pos ] =  sum(imfilter(center, mask))/(sub_img_size^2)
            mask_pos +=1
        end

        position += 1;
    end

    return local_correlation
end


"""
    extract_pixels_from_img(img, indicies_set, video_dim_tuple)

Takes every frame from @video_array and extracts pixels which indicies are in
@indicies_set, thus creating video with only chosen indicies.
"""
function extract_pixels_from_img(img, indicies_set, video_dim_tuple)
   rows = size(indicies_set,2)
   columns = size(indicies_set,2)
   video_length = video_dim_tuple[3]

   extracted_pixels = zeros(rows, columns, video_length)
   extracted_pixels[:,:,] =
                img[indicies_set[1,:],indicies_set[2,:]]

   return extracted_pixels
end

"""
Returns evenly distributed centers of size `image_size`
"""
function get_local_img_centers(points_per_dim, img_size, shift=0, sub_img_size=0 )
    # TODO Applied teproray solution here, so it works only for local gradients
    # start = 0
    # (points_per_dim>shift) ? start_ind = ceil(Int, points_per_dim/2)+ shift :
    #                         start=shift
    start_ind = ceil(Int, sub_img_size/2)
    min_va,  = findmin(img_size)
    last_ind = min_va - start_ind

    set = broadcast(floor, Int, range(start_ind, step=sub_img_size,  stop=last_ind))
    num_indexes = size(set,1)

    centers = Any[]
    for row = 1:num_indexes
        for column = 1:num_indexes
            push!(centers, CartesianIndex(set[row], set[column]))
        end
    end
    return centers
end


"""
    get_img_local_centers(img_size, sub_img_size=10)

Tiles the image of size `img_size` into square subimages of size `sub_img_size`
and returns vector with CartesianIndex coordinates of the subimages centre in
original image.

Takes smaller value from `img_size` and then divides it by `sub_img_size`.
Resulting value will be the number of returned subimages per dimension.
"""
function get_img_local_centers(img_size, sub_img_size=10)
    # TODO Applied teproray solution here, so it works only for local gradients
    min_va,  = findmin(img_size)
    if sub_img_size > min_va
        throw(BoundsError( "`sub_img_size` is bigger than original image."))
    end

    start_ind = ceil(Int, sub_img_size/2)
    last_ind = min_va - start_ind

    set = broadcast(floor, Int, range(start_ind, step=sub_img_size,  stop=last_ind))
    num_indexes = size(set,1)

    centers = Any[]
    for row = 1:num_indexes
        for column = 1:num_indexes
            push!(centers, CartesianIndex(set[row], set[column]))
        end
    end
    return centers
end


"""
    vectorize_img(video)

Rearrenges the video so that set of n frames (2D matrix varying in
time) the set of vectors is returned, in which each row is a pixel, and each
column is the value of the pixel in n-th frame.
"""
function vectorize_img(img)
    rows, columns = size(img)
    num_of_elements = rows*columns

    vectorized_img = zeros(num_of_elements)

    index = 1;
    for row=1:rows
        for column=1:columns
            vectorized_img[index] = img[row, column];
            index = index+1;
        end
    end

    return vectorized_img
end




"""
    get_video_mask(points_per_dim, video_dimensions; distribution="uniform", sorted=true, x=1, y=1)

Returns matrix of size @points_per_dim x 2 in which indicies of video frame are
stored. The indicies are chosen based one the @distribution argument. One option
is uniform distribution, the second is random distribution.

Uniform distribution: distance between the points in given dimension is the
 even, but vertical distance may be different from horizontal distance between points. This depends on the size of a frame in a image.

Random distribution: the distance between the points is not constant, because
the points are chosen randomly in the ranges 1:horizontal size of frame,
1:vertical size of frame. The returned values may be sorted in ascending order,
if @sorted=true.
"""
function get_video_mask(points_per_dim, video_dimensions;
                            distribution="uniform", sorted=true, patch_params)
    video_height, video_width,  = video_dimensions
    x=patch_params["x"]
    y=patch_params["y"]
    spread=patch_params["spread"]

    if x == 1
        x = Int64(floor(video_width/2))
        @warn "Given x is to close to the border. Seeting the value to " x
    elseif x < Int64(ceil(points_per_dim/2))
        x = Int64(ceil(points_per_dim/2))
        @warn "Given x is to close to the border. Seeting the value to " x
    elseif x > video_width-Int64(ceil(points_per_dim/2))
        x = video_width - Int64(ceil(points_per_dim/2))
        @warn "Given x is to close to the border. Seeting the value to " x
    end

    if y == 1
        y = Int64(floor(video_height/2))
        @warn "Given y is to close to the border. Seeting the value to " y
    elseif y < Int64(ceil(points_per_dim/2))
        y = Int64(ceil(points_per_dim/2))
        @warn "Given y is to close to the border. Seeting the value to " y
    elseif y > video_height-Int64(ceil(points_per_dim/2))
        y = video_height - Int64(ceil(points_per_dim/2))
        @warn "Given y is to close to the border. Seeting the value to " y
    end

    if spread*points_per_dim+x > video_width || spread*points_per_dim+y > video_height
        @warn "Given patch parameters might result in indicies exceeding frame size."
    end

    if distribution == "uniform"
        columns = points_per_dim
        rows = points_per_dim

        # +1 is used so that the number of points returned is as requested
        row_step = Int64(floor(video_height/rows))
        column_step = Int64(floor(video_width/columns))

        (video_height/row_step != points_per_dim) ? row_step+=1 : row_step
        (video_width/column_step !=
                                points_per_dim) ? column_step+=1 : video_width

        vertical_indicies = collect(1:row_step:video_height)
        horizontal_indicies = collect(1:column_step:video_width)

        vertical_indicies = reshape(vertical_indicies, (1,points_per_dim))
        horizontal_indicies = reshape(horizontal_indicies, (1,points_per_dim))

        indicies_set = [vertical_indicies; horizontal_indicies]
    elseif distribution == "random"
        vertical_indicies = rand(1:video_height,1, points_per_dim)
        horizontal_indicies = rand(1:video_width,1, points_per_dim)

        if sorted
            vertical_indicies = sort(vertical_indicies[1,:])
            horizontal_indicies = sort(horizontal_indicies[1,:])

            vertical_indicies = reshape(vertical_indicies, (1,points_per_dim))
            horizontal_indicies =
                              reshape(horizontal_indicies, (1,points_per_dim))
        end
        indicies_set = [vertical_indicies; horizontal_indicies]
    elseif distribution == "patch"
        indicies_set = [collect(1:spread:(spread*points_per_dim)).+x collect(1:spread:(spread*points_per_dim)).+y]'
    end

   return indicies_set
end


"""
    get_gabor_mask_set(;filt_size=25, σ=[2], angle_1=[0], λ=[15], γ=[0.2],
                            angle_2=[0], re_part=true, im_part=false)

Returns set of gabor filters generated with given parameters. Parameters are described
below. Function uses Kernel.gabor() from ImageFiltering.

# Arguments
- `filt_size=30` : controls the patch in which filter is created, not wavelet itself
- `σ=2` : controls the width of the waves and thus number of cycles per unit
- `angle_1=0` : is the rotation in (half)radians, θ = pi*(angle_1/180)
- `λ=15` : controls the number of waves within the window- higher values- less waves
- `γ=0.2` : is the aspect ratio; small values give long filters
- `angle_2=0` : phase, pi*(angle_2/180)
- `re_part::Bool`: determines if real part of the Gabor filter is returned; real
    part is normalized to be in range [0,1]
- `im_part::Bool`: determines if imaginary part of the Gabor filter is returned
    imaginary part is normalized to be in range [0,1]

if both `re_part` and `im_part` are true, then absolute value of complex number
    of form `re_part + im_part im` is returned (it is also normalized to range
    [0,1]).
"""
function get_gabor_mask_set(;filt_size=25, σ=[2], angle_1=[0], λ=[15], γ=[0.2],
                            angle_2=[0], re_part=true, im_part=false)

    kernels = Any[]
    for sigma = σ
        for angle1 = angle_1
            θ = pi*(angle1/180)
            for lambda in λ
                for gamma in γ
                    for angle2 in angle_2
                        ψ = pi*(angle2/180)
                        kernel = Kernel.gabor(filt_size, filt_size,
                                        sigma,
                                        θ,
                                        lambda,
                                        gamma,
                                        ψ)
                        if re_part && !im_part
                            kernel[1] .+= abs(findmin(kernel[1])[1])
                            kernel[1] ./= findmax(kernel[1])[1]
                            push!(kernels,Gray.((kernel[1])))

                        elseif im_part && !re_part
                            kernel[2] .+= abs(findmin(kernel[2])[1])
                            kernel[2] ./= findmax(kernel[2])[1]
                            push!(kernels,Gray.((kernel[2])))

                        else
                            @debug "Using abs(re(A)+im(A))"
                            result = abs.(kernel[1] + kernel[2]im);
                            result .+= abs(findmin(result)[1])
                            result ./= findmax(result)[1]
                            push!(kernels,Gray.())
                        end
                    end # angle2
                end # gamma
            end # lambda
        end # angle1
    end # sigmas
    return kernels
end



"""
    rearrange_filters_arr(im_filter; showing_number=-1)

Creates image with elements stored in `im_filters`. `showing_number` determines
how many of the element from `im_fiter` are displayed.

`im_filters` is an array with elements of type Matrix{Gray}.
"""
function rearrange_filters_arr(im_filter; showing_number=-1)
    mask_size = size(masks[1],1)
    masks_num = length(masks)
    if showing_number == -1 || showing_number > masks_num
        max_indeks = masks_num
    else
        max_indeks = showing_number
    end
    columns = Int(ceil(sqrt(masks_num)))
    rows= Int(ceil(masks_num/columns))

    all_filters = zeros(Gray, rows*mask_size, columns*mask_size)

    mask_index = 1
    for row in 1:rows
        start_row = (row-1)*mask_size+1
        row_range = start_row:(start_row+mask_size-1)

        for col = 1:columns
            start_col = (col-1)*mask_size+1
            col_range = start_col:(start_col+mask_size-1)
            if mask_index > max_indeks
                break
            else
                all_filters[row_range, col_range] = masks[mask_index]
                mask_index += 1
            end
        end
        if mask_index > max_indeks
            break
        end
    end
    return all_filters
end

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
    get_local_img_gradients(img, centers, sub_img_size)

Computes the gradients in the subimage, takes the mean of sum of absolute
values of both hotizontal and vertical gradients as a representative of a
subimage.
"""
function get_local_img_gradients(img, centers, sub_img_size)
    @debug "Entering get_local_gradients"
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size
    h, w, len = get_video_dimension(video_array)
    extracted_pixels = zeros(sub_img_size, sub_img_size)

    img_grad = imgradients(img, KernelFactors.ando3, "replicate")
    img_grad_abs = map(abs, img_grad[1]) + map(abs, img_grad[2])
    for index_x = 1:size(centers,2)
        c_x = centers[2, index_x]
        for index_y = 1:size(centers,2)
            c_y = centers[1, index_y]
            sub_img  = img_grad_abs[(c_x-half_size):(c_x+half_size),
                                    (c_y-half_size):(c_y+half_size)]

            extracted_pixels[index_x, index_y] =mean(sub_img)
        end
    end

    return extracted_pixels
end



"""
    get_local_img_correlations(img, centers, sub_img_size, shift)

Computes the correlation between the subimages and subimages shifted by values
from range -`shift`:`shift` and returns array of size
length(`centers`) x length(`centers`).

Each of the subimage is center around values stored in  @centers
"""
function get_local_img_correlations(img, centers, sub_img_size, shift)
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size + shift
    h, w = size(img)
    extracted_pixels = zeros(sub_img_size, sub_img_size)

    for index_x = 1:size(centers,2)
        c_x = centers[2, index_x]
        for index_y = 1:size(centers,2)
            c_y = centers[1, index_y]
            subimage = img[(c_x-half_range):(c_x+half_range),
                            (c_y-half_range):(c_y+half_range)]
            center = img[(c_x-half_size):(c_x+half_size), (c_y-half_size):(c_y+half_size)]

            for left_boundary = 1:(2*shift+1)
                for lower_boundary = 1:(2*shift+1)
                    corelation = center .* subimage[left_boundary:left_boundary+sub_img_size-1, lower_boundary:lower_boundary+sub_img_size-1]
                    corelation = sum(corelation)
                    extracted_pixels[index_x, index_y] += corelation
                end
            end
            extracted_pixels[index_x, index_y] /= 256*(sub_img_size^2)*(shift*2)^2
        end
    end

    return extracted_pixels
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

# import Makie
import VideoIO
 using StatsBase
 using Images
 using ImageFeatures
  # using TestImages
 using ImageDraw
 using CoordinateTransformations
 # using Makie
 using Logging

export get_video_array_from_file,
        get_video_dimension,
        get_video_mask,
        extract_pixels_from_video,
        extract_pixels_from_video,
        vectorize_video,
        get_pairwise_correlation_matrix,
        get_average_from_tiles,
        rotate_img_around_center,
        rotate_vid_around_center,
        export_images_to_vid,
        rotate_and_save_video,
        get_local_correlations,
        get_local_centers,
        get_local_gradients,
        normalize_to_01,
        shift_to_non_negative,
        plotimg;

"""
 get_video_array_from_file(video_name)

Returns array to which video frames are copied. Frames are in grayscale.

Function opens stream, then loads the file and gets all the frames from a
video.
 """
function get_video_array_from_file(video_name)
   video_streamer = VideoIO.open(video_name) # candle
   video_array = Vector{Array{UInt8}}(undef,0);
   video_file = VideoIO.openvideo(video_streamer, target_format=VideoIO.AV_PIX_FMT_GRAY8)

   while !eof(video_file)
      push!(video_array,reinterpret(UInt8, read(video_file)))
   end
   close(video_file)

   return video_array
end

"""
    get_video_dimension(video_array)

Returns the tuple which contains width, height and the number of the frames of
array in whcih video was loaded.
"""
function get_video_dimension(video_array)
   v_hei = size(video_array[1],1)
   v_wid = size(video_array[1],2)
   v_len = size(video_array,1)

   video_dim_tuple = (video_height=v_hei, video_width=v_wid, video_length=v_len)

   return video_dim_tuple
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
    extract_pixels_from_video(video_array, indicies_set, video_dim_tuple)

Takes every frame from @video_array and extracts pixels which indicies are in
@indicies_set, thus creating video with only chosen indicies.
"""
function extract_pixels_from_video(video_array, indicies_set, video_dim_tuple)
   rows = size(indicies_set,2)
   columns = size(indicies_set,2)
   video_length = video_dim_tuple[3]

   extracted_pixels = zeros(rows, columns, video_length)
   for frame_number in 1:video_length
      extracted_pixels[:,:,frame_number] =
                video_array[frame_number][indicies_set[1,:],indicies_set[2,:]]
   end

   return extracted_pixels
end


"""
    vectorize_video(video)

Rearrenges the video so that set of n frames (2D matrix varying in
time) the set of vectors is returned, in which each row is a pixel, and each
column is the value of the pixel in n-th frame.
"""
function vectorize_video(video)
    video_length = size(video, 3)
    rows = size(video,1)
    columns = size(video,2)

    number_of_vectors = rows*columns

    vectorized_video = zeros(number_of_vectors, video_length)

    index = 1;
    for row=1:rows
        for column=1:columns
            vectorized_video[index,:] = video[row, column,:];
            index = index+1;
        end
    end

    return vectorized_video
end



"""
    get_average_from_tiles(extracted_pixels_matrix, N)

Fnction takes a 3D array in which video is stored and splits every frame into
non overlaping tiles of size NxN. If size of @extracted_pixels_matrix is not
square of N, then only N^2 x N^2 matrix will be used for averaging.
"""
function get_average_from_tiles(extracted_pixels_matrix, N)
    # N = size(extracted_pixels,1)
    num_frames = size(extracted_pixels_matrix,3)
    mask_matrix = ones(N, N)
    result_matrix = zeros(N, N, num_frames)
    col_index = 1
    row_index = 1

    for frame = 1:num_frames
        for col = 1:N:N^2
            for row = 1:N:N^2
                result_matrix[mod(col,N), mod(row,N), frame] =
                        dot(extracted_pixels_matrix[col:(col+N-1),
                            row:(row+N-1), frame], mask_matrix) ./N^2
                row_index += 1
            end
            col_index += 1
        end
    end
    return result_matrix
end


"""
    rotate_vid_around_center(img, rotation = 5pi/6)

Function rotates a video around the center of the image by @rotation radians and
 the outout into matrix.
"""
function rotate_vid_around_center(src_vid_path,src_vid_name; rotation = 5pi/6)
    video_array = []
    video_src_strm = VideoIO.open(src_vid_path*src_vid_name)
    video_src = VideoIO.openvideo(video_src_strm,
                                        target_format=VideoIO.AV_PIX_FMT_GRAY8)

    while !eof(video_src)
      img = read(video_src)
      img = rotate_img_around_center(img, rotation)

      push!(video_array,img)
    end
    close(video_src)

  return video_array
end


"""
    export_images_to_vid(video_array, dest_file)

Exports set of images stored in @video_array to the dest_file.

"""
function export_images_to_vid(video_array, dest_file)
    @debug "Exporting set of images to file"
    fname = dest_file

    video_dimensions = get_video_dimension(video_array)
    h = video_dimensions.video_height
    w = video_dimensions.video_width
    nframes = video_dimensions.video_length
    overwrite=true
    fps=30
    options = ``
    ow = overwrite ? `-y` : `-n`

    open(`ffmpeg
            -loglevel warning
            $ow
            -f rawvideo
            -pix_fmt rgb24
            -s:v $(h)x$(w)
            -r $fps
            -i pipe:0
            $options
            -vf "transpose=0"
            -pix_fmt yuv420p
            $fname`, "w") do out
        for i = 1:nframes
            write(out, convert.(RGB{N0f8}, clamp01.(video_array[i])))
        end
    end
    @debug "Video was saved"
end


"""
    rotate_and_save_video(src_vid_path, src_vid_name, dest_vid_name;
                                                                rotation=5pi/6)

Fuction opens the @src_vid_name file, collects all the frames and then rotates
the frame aroung the center and saves new video as @dest_vid_name at
@src_vid_path.

Function was tested for following extensions;
    .mov

A solution for writing to a video file was taken from:
https://discourse.julialang.org/t/creating-a-video-from-a-stack-of-images/646/7
"""
function rotate_and_save_video(src_vid_path, src_vid_name, dest_vid_name, rotation=5pi/6)
    @debug src_vid_path src_vid_name dest_vid_name

    if !isfile(src_vid_path*src_vid_name)
        @warn "Source file at given path does not exist. Please give another name."
        return
    elseif isfile(src_vid_path*dest_vid_name)
        @warn "File with destination video name at src_video_path already exists. Please give another name."
        return
    end

    video_array = rotate_vid_around_ceter(src_vid_path, src_vid_name)
    @debug "Video was rotated"

    export_images_to_exist_vid(video_array, src_vid_path*dest_vid_name)
    @info "The file was created:\n  $fname"
end


function get_local_centers(points_per_dim, video_dimensions, shift=0, sub_img_size=0 )
    /# TODO Applied teproray solution here, so it works only for local gradients
    # start = 0
    # (points_per_dim>shift) ? start_ind = ceil(Int, points_per_dim/2)+ shift :
    #                         start=shift
    start_ind = ceil(Int, sub_img_size/2)
    min_va,  = findmin(video_dimensions)
    last_ind = min_va - start_ind

    set = broadcast(floor, Int, range(start_ind, stop=last_ind,  length=points_per_dim))
    centers = [set set]'
    return centers
end

import numpy as np 
from PIL import Image

# Define pixel dimensions in microns for the 3D stack
# x and y dimensions: 0.2075665 microns per pixel
# z dimension: 0.4870975 microns per pixel

# Initialize lists to store results
result = []
result2 = []

# Loop through each frame (0 to 53) in the time series
for Frame in range(0, 54):
    print('Processing frame:', Frame)
    
    # Load the 2D segmented and tracked cells and the corresponding 'z-stack projected and summed' image 
    tracked_cells = Image.open(f'tracked_cells/tracked_cells_resized_t{str(Frame).zfill(2)}.tif')
    summed_zstack = Image.open(f'Folder_with_3D_images_of_whole_tissue/movie_t{str(Frame).zfill(2)}.tif')
    
    # Convert the images to NumPy arrays
    Zstack_array = np.array(summed_zstack)
    tracked_cells_array = np.array(tracked_cells)
    
    # Initialize arrays to hold volume and count data for each cell
    volume = np.zeros(256**3)
    Count = np.zeros(256**3)
    
    # Initialize lists to track cell IDs and cells to be removed
    cell_id = []
    cell_to_remove = []
    
    # Iterate through each pixel in the tracked cells array
    for i in range(np.shape(tracked_cells_array)[0]):
        for j in range(np.shape(tracked_cells_array)[1]):

            # Calculate the unique cell ID based on RGB values
            N = (tracked_cells_array[i, j][2] * 256**2) + 
            (tracked_cells_array[i, j][1] * 256) + 
            (tracked_cells_array[i, j][0])
            
            # Check if the value in the z-stack image is below the threshold (1020.0)
            if Zstack_array[i, j] < 1020.0:
                if N not in cell_to_remove:
                    cell_to_remove.append(N)
            else:
                # Accumulate the volume for each cell ID if it is not to be removed
                if N not in cell_to_remove:
                    volume[N] += Zstack_array[i, j]  # Sum height values for each cell
                    Count[N] += 1
                    if N not in cell_id:
                        cell_id.append(N)
    
    # Calculate the cell volumes and store the results
    for k in cell_id:
        # frame number, cell ID, average height, volume in microns^3
        result.append([
            Frame, k, 
            (volume[k] * 0.4870975 / (Count[k] * 255)),  # Average height of the cell
            (volume[k] * 0.4870975 * 0.2075665 * 0.2075665 / 255)  # Volume of the cell in microns^3
        ])
    
    # Store cell IDs to be removed
    for x in cell_to_remove:
        result2.append([Frame, x])

# Save the results to a CSV file
np.savetxt('result_3D_corrected.csv', result, delimiter=',', fmt='%s')	

# Uncomment the line below if you want to save the removed cells data
# np.savetxt('cells_removed.csv', result2, delimiter=',', fmt='%s')

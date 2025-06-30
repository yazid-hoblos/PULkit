from PIL import Image
import os

# Paths to folders
folder_top = "panorama_cazy_plots"
folder_bottom = "plots"

# Select first 3 images from panorama_cazy_plots (sorted)
top_images_files = sorted([f for f in os.listdir(folder_top) if f.endswith('.png')])[:3]
# Select first 3 images from plots (sorted)
bottom_images_files = sorted([f for f in os.listdir(folder_bottom) if f.endswith('.png')])[:3]

# Load images
top_images = [Image.open(os.path.join(folder_top, f)) for f in top_images_files]
bottom_images = [Image.open(os.path.join(folder_bottom, f)) for f in bottom_images_files]

# Assume all images have the same size within their group
# But top and bottom group images can differ in size, so resize bottom images to top images width if needed
top_width, top_height = top_images[0].size
bottom_width, bottom_height = bottom_images[0].size

# Optional: Resize bottom images to match width of top images for neat alignment
if bottom_width != top_width or bottom_height != top_height:
    bottom_images = [img.resize((top_width, top_height), Image.ANTIALIAS) for img in bottom_images]

# Calculate combined image size (3 images per row, 2 rows)
combined_width = top_width * 3
combined_height = top_height * 2

# Create blank canvas
combined_img = Image.new('RGBA', (combined_width, combined_height))

# Paste top row images
for i, img in enumerate(top_images):
    combined_img.paste(img, (i * top_width, 0))

# Paste bottom row images
for i, img in enumerate(bottom_images):
    combined_img.paste(img, (i * top_width, top_height))

# Save and show
combined_img.save("report/combined_panorama_dbcan_vs_cazy_plots.png")
combined_img.show()

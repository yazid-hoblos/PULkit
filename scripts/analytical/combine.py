#!/usr/bin/env python3
"""
Image Grid Combiner

This script combines PNG images from two directories into a grid layout.
By default, it creates a 3x2 grid (3 images per row, 2 rows) using the first
3 images from each directory (sorted alphabetically).
"""

import os
import argparse
from PIL import Image
from pathlib import Path


def get_image_files(directory, count=3, extensions=('.png', '.jpg', '.jpeg')):
    """Get the first N image files from a directory (sorted)."""
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory '{directory}' not found")
    
    image_files = []
    for ext in extensions:
        image_files.extend([f for f in os.listdir(directory) if f.lower().endswith(ext)])
    
    if not image_files:
        raise ValueError(f"No image files found in '{directory}'")
    
    return sorted(image_files)[:count]


def load_and_resize_images(directory, image_files, target_size=None):
    """Load images and optionally resize them to target size."""
    images = []
    for filename in image_files:
        img_path = os.path.join(directory, filename)
        img = Image.open(img_path)
        
        if target_size:
            img = img.resize(target_size, Image.Resampling.LANCZOS)
        
        images.append(img)
    
    return images


def create_image_grid(top_images, bottom_images, output_path):
    """Create a grid layout with top and bottom row images."""
    if len(top_images) != len(bottom_images):
        raise ValueError("Number of top and bottom images must be equal")
    
    cols = len(top_images)
    rows = 2
    
    # Get dimensions from first image
    img_width, img_height = top_images[0].size
    
    # Calculate combined dimensions
    combined_width = img_width * cols
    combined_height = img_height * rows
    
    # Create blank canvas
    combined_img = Image.new('RGB', (combined_width, combined_height), color='white')
    
    # Paste top row images
    for i, img in enumerate(top_images):
        x_pos = i * img_width
        combined_img.paste(img, (x_pos, 0))
    
    # Paste bottom row images
    for i, img in enumerate(bottom_images):
        x_pos = i * img_width
        combined_img.paste(img, (x_pos, img_height))
    
    # Save the combined image
    combined_img.save(output_path)
    print(f"Combined image saved to: {output_path}")
    
    return combined_img


def main():
    """Main function to combine images from two directories."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("top_directory",
                       help="Directory containing images for the top row")
    parser.add_argument("bottom_directory",
                       help="Directory containing images for the bottom row")
    parser.add_argument("-o", "--output", default="report/combined_plots.png",
                       help="Output path for combined image (default: report/combined_plots.png)")
    parser.add_argument("-n", "--num-images", type=int, default=3,
                       help="Number of images to use from each directory (default: 3)")
    parser.add_argument("--resize", nargs=2, type=int, metavar=("WIDTH", "HEIGHT"),
                       help="Resize all images to specified width and height")
    parser.add_argument("--match-size", action="store_true",
                       help="Resize bottom images to match top image dimensions")
    parser.add_argument("--show", action="store_true",
                       help="Display the combined image after creation")
    
    args = parser.parse_args()
    
    try:
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(args.output)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Get image files from both directories
        print(f"Loading {args.num_images} images from '{args.top_directory}'...")
        top_files = get_image_files(args.top_directory, args.num_images, args.extensions)
        print(f"Top images: {top_files}")
        
        print(f"Loading {args.num_images} images from '{args.bottom_directory}'...")
        bottom_files = get_image_files(args.bottom_directory, args.num_images, args.extensions)
        print(f"Bottom images: {bottom_files}")
        
        # Load images
        top_images = load_and_resize_images(args.top_directory, top_files)
        bottom_images = load_and_resize_images(args.bottom_directory, bottom_files)
        
        # Determine target size for resizing
        target_size = None
        if args.resize:
            target_size = tuple(args.resize)
            print(f"Resizing all images to: {target_size}")
        elif args.match_size:
            target_size = top_images[0].size
            print(f"Resizing bottom images to match top size: {target_size}")
        
        # Resize if needed
        if target_size:
            if args.resize or args.match_size:
                top_images = load_and_resize_images(args.top_directory, top_files, target_size)
            if target_size:
                bottom_images = load_and_resize_images(args.bottom_directory, bottom_files, target_size)
        
        # Create combined image
        print("Creating combined image...")
        combined_img = create_image_grid(top_images, bottom_images, args.output)
        
        # Show image if requested
        if args.show:
            combined_img.show()
        
        print("Image combination complete!")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
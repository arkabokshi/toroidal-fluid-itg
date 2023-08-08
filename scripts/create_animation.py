"""
Creates an animated gif using several images as individual frames.

Requires the imageio library. See the docs here:
https://imageio.readthedocs.io/en/stable/index.html
"""
from pathlib import Path
import os
import imageio.v2 as imageio

images = []
for filename in sorted(Path("./images").iterdir(), key=os.path.getmtime):
    print("Reading image data from:", filename)
    images.append(imageio.imread(filename))
print("Saving .gif, please wait...")
imageio.mimsave("GPO_3D.gif", images, duration=125)  # duration is in ms (per frame)

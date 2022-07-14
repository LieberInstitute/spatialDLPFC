#   Image segmentation is extremely slow with squidpy on our full-resolution
#   images, and appears highly sensitive to environment variables such as
#   OMP_NUM_THREADS, but it isn't clear what values for these variables results
#   in best speeds. This script times segmentation on a test image.

import squidpy as sq
import numpy as np
import time

NUM_REPS = 10

start = time.perf_counter()

#   Load, preprocess and segment test image
for _ in range(NUM_REPS):
    img = sq.datasets.visium_fluo_image_crop()
    sq.im.process(img=img, layer="image", method="smooth")
    sq.im.segment(
        img=img,
        layer="image_smooth",
        method="watershed",
        channel=0,
    )

end = round(time.perf_counter() - start, 1)

print('Time elapsed for {} reps: {}s.'.format(NUM_REPS, end))

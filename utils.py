import colorsys
from PIL import Image

def format_image(image, expected_width):
    # try each rotation to find the finder pattern
    for angle in (0, 90, 180, 270):
        cand = image.rotate(angle, expand=True)
        # convert to grayscale for bar detection
        gray = cand.convert('L', dither=Image.NONE)
        w, h = gray.size
        px = gray.load()
        mid_y = h // 2
        # find first dark pixel on middle row
        start = next((x for x in range(w) if px[x, mid_y] < 128), None)
        if start is None:
            continue
        # find where the dark run ends
        end = next((x for x in range(start+1, w) if px[x, mid_y] >= 128), None)
        if end is None:
            continue
        bar_w = end - start
        # skip if there is not enough space for finder bars
        if start + 2*bar_w + bar_w//2 >= w:
            continue
        rgb = cand.convert('RGB').load()
        # sample the three finder bar centers
        centers = [
            rgb[start + bar_w//2,        mid_y],
            rgb[start + bar_w + bar_w//2, mid_y],
            rgb[start + 2*bar_w + bar_w//2, mid_y],
        ]
        def is_black(c): return sum(c)/3 < 64
        def is_colored(c):
            h_, s_, v_ = colorsys.rgb_to_hsv(c[0]/255, c[1]/255, c[2]/255)
            return s_ > 0.5 and v_ > 0.5
        # confirm pattern black colored black
        if is_black(centers[0]) and is_colored(centers[1]) and is_black(centers[2]):
            rotation = angle
            detected_bar_width = bar_w
            break
    else:
        raise ValueError("Could not detect rotation/finder pattern")

    # rotate image to correct orientation
    image = image.rotate(rotation, expand=True, resample=Image.Resampling.NEAREST)
    width, height = image.size
    # choose bar width for cropping

    rescale_ratio = expected_width / width

    # rescale image to expected width and height
    image = image.resize((round(width * rescale_ratio), round(height * rescale_ratio)), Image.Resampling.NEAREST)
    


    # save for inspection
    image.save('barcode_formatted.png')
    return image


def get_array_from_image(image, barWidth, barHeight):
    # convert image to RGB for sampling
    image = image.convert('RGB')
    w, h = image.size
    n = w // barWidth
    px = image.load()
    array = bytearray(n)

    # sample each bar by averaging all pixels in its column block
    for i in range(n):
        r_sum = g_sum = b_sum = 0
        x0 = i * barWidth
        # iterate every pixel in this bar’s width and full height
        for dx in range(barWidth):
            for y in range(h):
                r, g, b = px[x0 + dx, y]
                r_sum += r
                g_sum += g
                b_sum += b
        count = barWidth * h
        r_avg = r_sum // count
        g_avg = g_sum // count
        b_avg = b_sum // count

        # convert average color to HSV
        h_, s_, v_ = colorsys.rgb_to_hsv(r_avg/255.0,
                                         g_avg/255.0,
                                         b_avg/255.0)
        # decide if bar is black or colored
        if v_ < 0.1 or (v_ > 0.9 and s_ < 0.1):
            array[i] = 0
        else:
            array[i] = int(round(h_ * 255))

    # remove finder and logo gap as before
    logo_width_in_bars = barHeight // barWidth
    return array[:3] + array[3 + logo_width_in_bars:]

    # we know 3 bars are finder, and then barHeight/barWidth bars are logo gap
    # remove logo gap from array
    logo_width_in_bars = barHeight // barWidth
    array_no_logo = array[:3] + array[3 + logo_width_in_bars:]
    return array_no_logo

def get_data_from_array(array, barWidth):
    # read error correction level from second bar
    error_correction_level = array[1]
    # strip finder bars and logo gap and end bar
    data = array[3  : -1]
    idx = data.find(b'\x00\x00')
    if idx != -1:
        data = data[:idx]
    return data, error_correction_level

def compare_results(original, decoded):
    # count mismatched characters
    wrong = sum(1 for i in range(len(original))
                if i >= len(decoded)
                or original[i] != decoded[i])
    return wrong

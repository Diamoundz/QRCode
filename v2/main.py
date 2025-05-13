from PIL import Image, ImageDraw
import random
import colorsys
import reedsolo
from reedsolo import RSCodec, ReedSolomonError

# Barcode data includes finder pattern: 3 start + N data + 1 end
arraySize = 1000
barWidth = 10              # width of each bar in pixels
error_correction_level = 80  # Example error correction level (0–255)

def set_finder_pattern(array, error_correction_level=0):
    array[0] = 0                       # black
    array[1] = error_correction_level # colored (encodes ECC level)
    array[2] = 0                       # black
    array[-1] = 0                      # black

def get_array_image(array, barWidth):
    """
    Draw the barcode.
    Finder bars (0, 2, last) are always black.
    ECC bar (index 1) is always colored.
    All OTHER positions:
        byte == 0   → white
        byte == 255 → white
        else        → colored by hue
    """
    img_width = len(array) * barWidth
    img_height = 100
    img = Image.new('RGB', (img_width, img_height), 'white')
    draw = ImageDraw.Draw(img)

    n = len(array)
    for i, byte in enumerate(array):
        # --- finder pattern ---
        if i in (0, 2, n - 1):
            rgb = (0, 0, 0)
        # --- ECC bar ---
        elif i == 1:
            hue = (byte / 255.0) * 360
            rgb = tuple(int(c * 255) for c in
                        colorsys.hsv_to_rgb(hue / 360.0, 1.0, 1.0))
        # --- data area ---
        elif byte in (0, 255):
            rgb = (255, 255, 255)
        else:
            hue = (byte / 255.0) * 360
            rgb = tuple(int(c * 255) for c in
                        colorsys.hsv_to_rgb(hue / 360.0, 1.0, 1.0))

        x0 = i * barWidth
        draw.rectangle([x0, 0, x0 + barWidth, img_height], fill=rgb)

    return img

def save_image(image, filename):
    image.save(filename)

def load_image(filename):
    return Image.open(filename)

def format_image(image, targetBarWidth=None):
    """
    1) Try each 90° rotation (with expand=True) to find the 0-1-0 finder.
    2) Once we know the correct angle & bar_w, center-crop back to the original size.
    3) Finally, crop to an integer number of bars.
    """
    import colorsys

    orig_w, orig_h = image.size

    # --- Rotation detection ---
    for angle in (0, 90, 180, 270):
        cand = image.rotate(angle, expand=True)
        gray = cand.convert('L')
        w, h = gray.size
        px = gray.load()

        mid_y = h // 2
        start = next((x for x in range(w) if px[x, mid_y] < 128), None)
        if start is None:
            continue
        end = next((x for x in range(start+1, w) if px[x, mid_y] >= 128), None)
        if end is None:
            continue
        bar_w = end - start

        # ensure 3 finder bars fit
        if start + 2*bar_w + bar_w//2 >= w:
            continue

        rgb = cand.convert('RGB').load()
        centers = [
            rgb[start + bar_w//2,        mid_y],
            rgb[start + bar_w + bar_w//2, mid_y],
            rgb[start + 2*bar_w + bar_w//2, mid_y],
        ]

        def is_black(c): return sum(c)/3 < 64
        def is_colored(c):
            h_, s_, v_ = colorsys.rgb_to_hsv(c[0]/255, c[1]/255, c[2]/255)
            return s_ > 0.5 and v_ > 0.5

        if is_black(centers[0]) and is_colored(centers[1]) and is_black(centers[2]):
            rotation = angle
            detected_bar_width = bar_w
            break
    else:
        raise ValueError("Could not detect rotation/finder pattern")

    # --- Recover a full-size, correctly rotated image ---
    cand = image.rotate(rotation, expand=True)
    # center-crop back to orig_w × orig_h
    left = (cand.width  - orig_w) // 2
    top  = (cand.height - orig_h) // 2
    image = cand.crop((left, top, left + orig_w, top + orig_h))
    width, height = orig_w, orig_h

    # --- Final cropping to whole bars ---
    bar_width = targetBarWidth or detected_bar_width
    total_bars = width // bar_width
    new_width  = total_bars * bar_width
    formatted   = image.crop((0, 0, new_width, height))

    return formatted, bar_width


def get_array_from_image(image, barWidth):
    """
    Decode a barcode by sampling a single horizontal line (mid-height)
    and averaging all pixels across each bar’s width.
    """
    image = image.convert('RGB')
    w, h = image.size
    n = w // barWidth
    px = image.load()
    y = h // 2

    array = bytearray(n)
    for i in range(n):
        r_sum = g_sum = b_sum = 0
        x0 = i * barWidth
        for dx in range(barWidth):
            r, g, b = px[x0 + dx, y]
            r_sum += r; g_sum += g; b_sum += b
        r_avg = r_sum // barWidth
        g_avg = g_sum // barWidth
        b_avg = b_sum // barWidth

        # convert to HSV
        h_, s_, v_ = colorsys.rgb_to_hsv(r_avg/255.0,
                                         g_avg/255.0,
                                         b_avg/255.0)
        if v_ < 0.1 or (v_ > 0.9 and s_ < 0.1):
            array[i] = 0
        else:
            array[i] = int(round(h_ * 255))
    return array

def get_data_from_array(array):
    error_correction_level = array[1]
    data = array[3:-1]
    return data, error_correction_level

def simulate_image_noise(image, noise_level=0.1):
    """
    Add pixel noise and random 90°-step rotation.
    """
    pixels = image.load()
    width, height = image.size
    for _ in range(int(width * height * noise_level)):
        x = random.randint(0, width - 1)
        y = random.randint(0, height - 1)
        pixels[x, y] = (random.randint(0, 255),
                        random.randint(0, 255),
                        random.randint(0, 255))
    rotation = random.choice([0, 90, 180, 270])
    print("Rotation:", rotation)
    return image.rotate(rotation, expand=True)

def compare_results(original, decoded):
    wrong = sum(1 for i in range(len(original)) if original[i] != decoded[i])
    return wrong

def max_encodable_message_bytes(bits, nsym):
    capacity_bytes = bits // 8
    return max(0, capacity_bytes - nsym*2)

def encode_data(data: bytes, array_size: int, error_correction_level: int) -> bytearray:
    rsc = RSCodec(error_correction_level)
    capacity = array_size - 4
    max_payload = capacity - rsc.nsym
    if len(data) > max_payload:
        data = data[:max_payload]
    return bytearray(rsc.encode(data))

def decode_data(codeword, error_correction_level):
    rsc = RSCodec(error_correction_level)
    msg = rsc.decode(codeword)
    if isinstance(msg, tuple):
        msg = msg[0]
    return bytearray(msg)

def run_test():
    """
    Wrap the existing Test Run logic in a try/except for cleaner error reporting.
    """
    try:
        # Create full array including finder pattern
        byteArray = bytearray(arraySize)
        data_str = "Hello, World! This is a rather long string - Once upon a time..."
        byteArray[3:-1] = data_str.encode('utf-8')
        print("Original data (str):", data_str)
        set_finder_pattern(byteArray, error_correction_level)

        encoded = encode_data(byteArray[3:-1], arraySize, error_correction_level)
        print("Encoded data length:", len(encoded))
        encoded_with_finder = bytearray(arraySize)
        encoded_with_finder[3:-1] = encoded
        set_finder_pattern(encoded_with_finder, error_correction_level)

        barcode_image = get_array_image(encoded_with_finder, barWidth)
        barcode_image = simulate_image_noise(barcode_image, noise_level=0.0)
        save_image(barcode_image, 'barcode.png')

        loaded_image = load_image('barcode.png')
        formatted_image, new_bar_width = format_image(loaded_image, targetBarWidth=barWidth)

        read_array = get_array_from_image(formatted_image, new_bar_width)
        data_chunk, ecc = get_data_from_array(read_array)
        print("decoded ecc level:", ecc)

        decoded = decode_data(data_chunk[:len(encoded)], ecc)
        wrong = compare_results(byteArray[3:-1][:len(decoded)], decoded)
        print("Decoded data (str):", decoded.decode('latin-1'))
        print("Number of wrong bytes:", wrong)

    except (ValueError, ReedSolomonError) as e:
        print("Error during test run:", e)

if __name__ == "__main__":
    run_test()

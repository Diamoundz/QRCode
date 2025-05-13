from PIL import Image, ImageDraw
import colorsys
import random
import math

def set_finder_pattern(array, error_correction_level=0):
    array[0] = 0
    array[1] = error_correction_level
    array[2] = 0
    array[-1] = 0

def get_array_image(array, barWidth, img_height):
    arrayCpy = array.copy()

    logo = Image.open('logo.jpg')
    # compute scaling factor to match barcode height
    w0, h0 = logo.size
    scale = img_height / min(w0, h0)
    #print logo size
    # resize while preserving aspect ratio
    logo = logo.resize((int(w0 * scale), int(h0 * scale)), Image.Resampling.LANCZOS)
    # crop to square centered on resized logo
    lw, lh = logo.size
    left = (lw - img_height) // 2
    top  = (lh - img_height) // 2
    logo = logo.crop((left, top, left + img_height, top + img_height))


    # compute how many bars the logo occupies
    logo_width_in_bars = max(1, round(logo.width / barWidth))
    # insert empty bars into data array after finder pattern, insert() arguments are (index, value)
    for _ in range(logo_width_in_bars):
        arrayCpy.insert(3, 0)

    # compute full image width based on updated array
    img_width = len(arrayCpy) * barWidth
    # create blank barcode image
    img = Image.new('RGB', (img_width, img_height), 'white')
    draw = ImageDraw.Draw(img)
    n = len(arrayCpy)

    # draw each bar according to array value
    for i, byte in enumerate(arrayCpy):
        if i in (0, 2, n - 1):
            rgb = (0, 0, 0)
        elif i == 1:
            hue = (byte / 255.0) * 360
            rgb = tuple(int(c * 255) for c in colorsys.hsv_to_rgb(hue / 360.0, 1.0, 1.0))
        elif byte in (0, 255):
            rgb = (255, 255, 255)
        else:
            hue = (byte / 255.0) * 360
            rgb = tuple(int(c * 255) for c in colorsys.hsv_to_rgb(hue / 360.0, 1.0, 1.0))
        x0 = i * barWidth
        draw.rectangle([x0, 0, x0 + barWidth, img_height], fill=rgb)

    # paste the prepared logo after the three finder bars
    x_position = barWidth * 3
    img.paste(logo, (x_position, 0))

    return img

def save_image(image, filename):
    image.save(filename)

def load_image(filename):
    return Image.open(filename)

def simulate_image_noise(image, noise_level=0.1):
    # get pixel map for random noise insertion
    pixels = image.load()
    width, height = image.size
    # corrupt a fraction of pixels with random colors
    for _ in range(int(width * height * noise_level)):
        x = random.randint(0, width - 1)
        y = random.randint(0, height - 1)
        pixels[x, y] = (random.randint(0, 255),
                        random.randint(0, 255),
                        random.randint(0, 255))
    # apply a random rotation for testing
    rotation = random.choice([0, 90, 180, 270])

    # Scale the image to a random size
    
    scale_factor = random.uniform(0.8, 1.2)
    #scale_factor = 1
    new_width = int(image.width * scale_factor)
    new_height = int(image.height * scale_factor)
    image = image.resize((new_width, new_height), Image.Resampling.NEAREST)
    return image.rotate(rotation, expand=True)

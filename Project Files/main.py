# import modules from barcode system
from barcode_drawer import set_finder_pattern, get_array_image, save_image, load_image, simulate_image_noise
from encoder import encode_data, decode_data
from utils import format_image, get_array_from_image, get_data_from_array, compare_results
from reedsolo import ReedSolomonError # type: ignore

barHeight = 100
arraySize = 150
barWidth = 20
error_correction_level = 80
noise_level = 0.01
input_string = "Hello, World! This is a rather long string - Once upon a time..."

def configure_test():
    global arraySize, barWidth, error_correction_level, input_string, noise_level
    print("Configure test parameters ===============")
    try:
        while True:
            # ask user for array size
            arraySize = int(input("Enter the size of the byte array (default 150): ") or 150)
            if arraySize <= 0:
                print("Array size must be a positive integer.")
                continue

            # ask user for bar width
            barWidth = int(input("Enter the width of the bars (default 20): ") or 20)
            if barWidth < 1:
                print("Bar width must be a positive integer.")
                continue

            # ask user for error correction level
            error_correction_level = int(input("Enter the error correction level (default 80): ") or 80)
            if error_correction_level < 0 or error_correction_level > 255:
                print("Error correction level must be between 0 and 255.")
                continue

            # ask user for noise level
            noise_level = float(input("Enter the noise level (default 0.01): ") or 0.01)
            if noise_level < 0.0 or noise_level > 1.0:
                print("Noise level must be between 0.0 and 1.0.")
                continue

            # ask user for data string to encode
            new_input = input("Enter the string to encode (default is preset text): ")
            if new_input:
                input_string = new_input
            if len(input_string) == 0:
                print("Input string cannot be empty.")
                continue

            break

    except ValueError as e:
        # use defaults if input is invalid
        print(f"Invalid input: {e}. Using default values.")

def run_test():
    try:
        # prepare raw byte array
        byteArray = bytearray(arraySize)
        data_bytes = input_string.encode('latin-1')

        payload_size = arraySize - 4
        # fill payload with input data

        print("Original data (str):", input_string)

        # insert finder pattern into array
        set_finder_pattern(byteArray, error_correction_level)
        # encode data with error correction
        encoded = encode_data(data_bytes, arraySize, error_correction_level)
        # prepare final array with finder pattern
        encoded_with_finder = bytearray(arraySize)
        encoded_with_finder[3:-1] = encoded

        set_finder_pattern(encoded_with_finder, error_correction_level)

        # generate barcode image from array
        barcode_image = get_array_image(encoded_with_finder, barWidth, barHeight)
        # apply noise to image for testing
        barcode_image = simulate_image_noise(barcode_image, noise_level)


        save_image(barcode_image, 'barcode.png')

        # ========================================== Now we simulate a test run where we read the image back =================================

        loaded_image = load_image('barcode.png')

        # detect and align barcode orientation
        formatted_image = format_image(loaded_image, expected_width= arraySize * barWidth  + barHeight)
        # convert image back to array
        read_array = get_array_from_image(formatted_image, barWidth, barHeight)
        # extract data and error correction level from array
        data_chunk, ecc = get_data_from_array(read_array, barWidth)

        # decode data with error correction
        decoded = decode_data(data_chunk[:len(encoded)], ecc)

        decoded_str = decoded.decode('latin-1').rstrip('\x00')
        wrong = compare_results(input_string, decoded_str)

        print("Decoded data (str):", decoded_str)
        print("Number of wrong bytes:", wrong)

    except (ValueError, ReedSolomonError) as e:
        print("Error during test run:", e)

if __name__ == "__main__":
    configure_test()
    run_test()

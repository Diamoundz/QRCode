from reedsolo import RSCodec  # type: ignore



def encode_data(data: bytes, array_size: int, error_correction_level: int) -> bytearray:
    # create codec with given ecc level
    rsc = RSCodec(error_correction_level)
    # number of ecc bytes produced
    actual_nsym = rsc.nsym
    capacity = array_size - 4

    # maximum payload after reserving ecc bytes
    max_payload = capacity - actual_nsym

    # ensure there is space for at least one byte of data
    if max_payload <= 0:
        raise ValueError(
            f"Byte array too small ({array_size}) for ECC level {error_correction_level} "
            f"(requires at least {actual_nsym + 1 + 4} total size)."
        )

    # truncate data if it exceeds available payload
    if len(data) > max_payload:
        print(f"Warning: Data truncated to {max_payload} bytes.")
        print("Consider increasing array_size, reducing error_correction_level, or reducing data size.")
        data = data[:max_payload]

    # encode payload and ecc bytes
    encoded = bytearray(rsc.encode(data))
    #pad the rest with 0x00 until the array size is reached
    if len(encoded) < array_size-4:
        encoded += bytearray(array_size - 4  - len(encoded))

    return encoded
    
    

# decode codeword bytes using reed solomon
def decode_data(codeword, error_correction_level):
    # create codec with same ecc level
    rsc = RSCodec(error_correction_level)
    # decode codeword and extract message
    msg = rsc.decode(codeword)

    if isinstance(msg, tuple):
        msg = msg[0]
    return bytearray(msg)

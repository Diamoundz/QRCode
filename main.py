from matrix import Matrix




def main():
    message = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etiam interdum dui porttitor dui cursus, at laoreet nulla gravida. Pellentesque ac Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etiafaucibus purus. Mauris mollis tincidunt lacus a bibendum. Integer Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. EtiavelLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdumLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia posuere. Curabitur sed sollicitudin diam. Etia facilisis turpis. Aliquam eleifend scelerisque nibh, ac elementum diam egestas vitae. Pellentesque id faucibus odio. Lorem ipsum dolor sit amet, consecteturLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. EtiaLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. EtiaLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia adipiscing elit. Phasellus euismod nunc nec purus hendrerit dictum. Integer nunc ligulaLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia, efficitur in lacus ut, ullamcorper pellentesque est. Vivamus eget gravida quam. AliquamLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia ut magna sit amet justo posuere maximus. Nunc ullamcorper nisi ac purus consequat, sagittis lobortis ligula pharetra. Suspendisse viverra sagittis nunc, rutrum tinciduntLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. EtiaLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia orci cursus sit amet. Interdum et malesuada fames ac anteLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia ipsum primis in faucibus. Curabitur placerat dui at mauris aliquam, et euismod odio porttitor. NamLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia erat tellus, ultricies vitae lectus ac, mollis blandit sem. Pellentesque dictum sapien a faucibus varius. Nam eget eratLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia nibh. Curabitur semper ante id diam bibendum tempus. Curabitur non fermentumLorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla vestibulum sagittis magna interdum posuere. Curabitur sed sollicitudin diam. Etia est. Sed nec lectus tempus velit varius consectetur eu ut ex. Cras mattis, dui eu lobortis rutrum, velit mauris elementum leo, vitae consectetur massa quam non ante. Mauris lacus erat, iaculis vel hendrerit eget, consequat nec nunc. Donec ullamcorper vel turpis faucibus commodo. In mollis iaculis nulla, a volutpat dolor scelerisque in."
    bits = encodeMessage(message)
    matrix = Matrix()
    truncatedLength = matrix.getMessageMaxSize()
    bits = bits[:truncatedLength]
    
    print("encoded bits:" + message[:truncatedLength//8])
    
    matrix = matrix.setBits(bits)
    
    encodedBits = matrix.values
    
    print("decoded bits:")
    print(matrix.decodeFromBitMatrix(encodedBits))
    
    
    matrix.displayMatrix()
    print()


def encodeMessage(message):
    # Encode 1 letter on a byte of 8 bits ASCII, return a list of bits for the message
    # Example: 'A' -> 01000001 (65 in decimal)
    # We return a string of bits containing the ASCII code of the message
    return ''.join(format(ord(char), '08b') for char in message)
    

main()
    
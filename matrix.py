import math
from PIL import Image, ImageDraw

def draw_hexagon(draw, center, side, fill_color):
    cx, cy = center
    angle_deg = 60
    angle_rad = math.pi / 180 * angle_deg
    vertices = [
        (
            cx + side * math.cos(angle_rad * i),
            cy + side * math.sin(angle_rad * i)
        )
        for i in range(6)
    ]
    draw.polygon(vertices, fill=fill_color)


class Matrix:
    SIZE = 25  # Defines a hexagon of radius (SIZE-1) around the center
    # this size needs to be > NSYM * 8 to have a non-null payload

    def __init__(self, finder_size=3):
        """
        Initialize the matrix with a dynamic finder pattern size.
        finder_size defines the offset from the hexagon boundary for the finder patterns.
        """
        self.finder_size = finder_size
        self.values = {}

    def _reserved_coords(self, reserved_levels=3):
        """
        Compute reserved coordinates for the finder patterns.
        TODO : Also add timing
        These coordinates are removed from the data area.
        """
        radius = self.SIZE - 1
        reserved = set()
        # Finder pattern centers (adjusted inward by 3 modules)
        finder_centers = [
            (-radius + self.finder_size, 0),            
            (radius - self.finder_size, -radius + self.finder_size),     
            (0, -radius + self.finder_size)             
        ]
        for center_q, center_r in finder_centers:
            for dq in range(-reserved_levels, reserved_levels + 1):
                for dr in range(max(-reserved_levels, -dq - reserved_levels),
                                min(reserved_levels, -dq + reserved_levels) + 1):
                    reserved.add((center_q + dq, center_r + dr))
        return reserved

    def getMessageMaxBits(self):
        # Total cells in a hexagon of radius (SIZE-1)
        radius = self.SIZE - 1
        total = 1 + 3 * radius * (radius + 1)
        # Subtract reserved finder pattern cells
        reserved = len(self._reserved_coords())
        return total - reserved  # available for data

    @classmethod
    def decodeFromBitMatrix(cls, bit_matrix):
            """
            Given a dictionary representing a hexagonal matrix of bits (with keys as axial coordinates (q, r)
            and values as 0 or 1), reconstitute the Matrix instance and decode the embedded message.
            
            The method deduces the hexagon's size from the keys and uses the default finder pattern size.
            Reserved finder pattern cells (computed from the deduced SIZE and finder_size) are skipped during decoding.
            """
            # Deduce the hexagon radius from the keys.
            radius = max(max(abs(q), abs(r), abs(-q - r)) for q, r in bit_matrix.keys())
            SIZE = radius + 1
            # Create a new instance with the default finder_size (3).
            instance = cls()
            instance.SIZE = SIZE
            # Populate instance.values from the provided bit_matrix (convert to booleans).
            instance.values = {coords: bool(bit_matrix[coords]) for coords in bit_matrix}
            
            # Read data cells in the same axial order as setBits, skipping reserved finder cells.
            bits = ""
            reserved_coords = instance._reserved_coords()
            for q in range(-radius, radius + 1):
                r1 = max(-radius, -q - radius)
                r2 = min(radius, -q + radius)
                for r in range(r1, r2 + 1):
                    if (q, r) in reserved_coords:
                        continue
                    bits += '1' if instance.values.get((q, r), False) else '0'
                    
            # Convert every 8 bits to a character.
            chars = []
            for i in range(0, len(bits), 8):
                byte = bits[i:i+8]
                if len(byte) == 8:
                    chars.append(chr(int(byte, 2)))
            final = ''.join(chars)
            return final.encode('latin-1')

    def setBits(self, bits):
        radius = self.SIZE - 1
        reserved_coords = self._reserved_coords()  # compute reserved positions

        idx = 0
        # Iterate over all axial coordinates in the overall hexagon
        for q in range(-radius, radius + 1):
            r1 = max(-radius, -q - radius)
            r2 = min(radius, -q + radius)
            for r in range(r1, r2 + 1):
                if (q, r) in reserved_coords:
                    # Mark reserved positions with None so they won't display data.
                    self.values[(q, r)] = None
                else:
                    if idx < len(bits):
                        self.values[(q, r)] = bits[idx] == '1'
                        idx += 1
                    else:
                        self.values[(q, r)] = False
        #print how many bits were set
        
        return self

    def draw_position_hexagon(self, draw, center_q, center_r, side, levels=3):
        # Draw finder pattern cells (both outer ring and center) as per our spec.
        for dq in range(-levels, levels + 1):
            for dr in range(max(-levels, -dq - levels),
                            min(levels, -dq + levels) + 1):
                # Only draw cells for the outer ring (dist==levels) and center (dist==0)
                dist = max(abs(dq), abs(dr), abs(-dq - dr))
                if dist == levels or dist == 0:
                    q = center_q + dq
                    r = center_r + dr
                    x = self.x_offset + 1.5 * side * q
                    y = self.y_offset + math.sqrt(3) * side * (r + q / 2)
                    # Finder pattern: outer ring and center are black
                    color = (0, 0, 0)
                    draw_hexagon(draw, (x, y), side, fill_color=color)

    def draw_timing_strip(self, draw, left_q, left_r, right_q, right_r):
        # Draw a timing strip between two points in the hexagonal grid.
        # The strip is a line of alternating colors (black and white).
        x1 = self.x_offset + 1.5 * left_q * 30
        y1 = self.y_offset + math.sqrt(3) * left_r * 30 + 15 * left_q
        x2 = self.x_offset + 1.5 * right_q * 30
        y2 = self.y_offset + math.sqrt(3) * right_r * 30 + 15 * right_q
        # Draw the line with alternating colors
        for i in range(0, int(math.hypot(x2 - x1, y2 - y1)), 10):
            color = (0, 0, 0) if i % 20 == 0 else (255, 255, 255)
            x = int(x1 + (x2 - x1) * i / math.hypot(x2 - x1, y2 - y1))
            y = int(y1 + (y2 - y1) * i / math.hypot(x2 - x1, y2 - y1))
            draw_hexagon(draw, (x, y), 15, fill_color=color)
            

    def displayMatrix(self, pixel_size=30, output_file='matrix_hexagonal.png'):
        side = pixel_size
        hex_height = math.sqrt(3) * side
        N = self.SIZE

        # Calculate overall image dimensions to fit the full hexagon
        img_width = int((2 * N - 1) * 1.5 * side + 2 * side)
        img_height = int((2 * N - 1) * hex_height + 2 * side)
        img = Image.new('RGB', (img_width, img_height), color='white')
        draw = ImageDraw.Draw(img)

        # Set offsets to center the hexagon.
        self.x_offset = img_width // 2
        self.y_offset = img_height // 2

        # Draw main hex grid for non-reserved cells.
        for (q, r), value in self.values.items():
            x = self.x_offset + 1.5 * side * q
            y = self.y_offset + hex_height * (r + q / 2)
            # If value is None (reserved for finder pattern), force white background.
            if value is None:
                color = (255, 255, 255)
            else:
                color = (0, 0, 0) if value else (255, 255, 255)
            draw_hexagon(draw, (x, y), side, fill_color=color)

        # Draw finder patterns so they stand out.
        radius = self.SIZE - 1
        self.draw_position_hexagon(draw, -radius + self.finder_size, 0, side)         # top-left finder
        self.draw_position_hexagon(draw, radius - self.finder_size, -radius + self.finder_size, side)  # top-right finder
        self.draw_position_hexagon(draw, 0, -radius + self.finder_size, side)    
        
        self.draw_timing_strip(draw, 0,0, -radius + self.finder_size, 0)  # horizontal timing strip
        
        img.save(output_file)
        img.show()

# Example usage:
if __name__ == '__main__':
    # Adjusted test bits: only fill available (non-reserved) data bits.
    matrix = Matrix()
    max_bits = matrix.getMessageMaxBits()
    # For testing, alternate ones and zeros.
    test_bits = ''.join(['1' if i % 2 == 0 else '0' for i in range(max_bits)])
    matrix.setBits(test_bits).displayMatrix(pixel_size=20)

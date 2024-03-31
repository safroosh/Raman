import re

def main():
    file_name = 'scf.out'
    siesta_cor_file = 'siesta.cor'

    try:
        with open(file_name, 'r') as file:
            content = file.read()

        num_of_species = int(re.search(r'NumberOfSpecies\s+(\d+)', content).group(1))
        num_of_atoms = int(re.search(r'NumberOfAtoms\s+(\d+)', content).group(1))
        print(num_of_atoms)

        lattice = extract_lattice_vectors(content)
        print("Lattice Vectors:", lattice)

        elements, counts = extract_elements_and_counts(content)
        print(elements)
        print(counts)
        
        target_line = "outcoor: Atomic coordinates (Ang)"
       
        
        found_lines = find_and_read_lines_from_string(content, target_line, num_of_atoms)
       # print(found_lines)
       
        coordinates_by_element = extract_coordinates_by_element(found_lines)
        
        with open("diagon.in", 'w') as file:
            file.write(f"{len(elements)} # number of species\n")
            file.write(f"{sum(counts)} # number of displaced atoms\n")
            file.write(" ".join(elements)+" #  chemical elements\n")
            file.write(" ".join(map(str, counts)) +" # number of atoms of each element\n")
            file.write("0.2 # scaling amplitude, 0.2 is good for Raman, change to 1.0 for visualization")

       
        with open(siesta_cor_file, 'w') as file:
           
    # Iterate through the dictionary
           for element, coordinates in coordinates_by_element.items():
               
        # Write the element symbol
               file.write(f"{element} {len(coordinates)}\n")
        # Iterate through the list of coordinates for this element
               for coords in coordinates:
                   
            # Format and write the coordinates
                   formatted_coords = '   '.join(coords)
                   file.write(f"{formatted_coords}\n")


    except FileNotFoundError:
        print(f"File {file_name} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def extract_lattice_vectors(content):
    # Dummy function for extracting lattice vectors
    
    lattice_vectors = re.findall(r'%block LatticeVectors\n([\d\.\s-]+)', ''.join(content), re.MULTILINE)
    lattice = [list(map(float, vec.split())) for vec in lattice_vectors[0].split('\n') if vec.strip() != '']
    return lattice

def extract_elements_and_counts(content):
    chem_spec_label_block = re.findall(r'%block Chemical_Species_Label\n([\d\.\s\w-]+)\n%endblock Chemical_Species_Label', ''.join(content), re.MULTILINE)
    chems = [line.strip() for line in chem_spec_label_block[0].split('\n') if line.strip() != '']

    elements = [line.split()[-1] for line in chems]
    natoms = [int(line.split()[1]) for line in chems]

    return (elements, natoms)


def extract_coordinates(lines):
    coordinates = []
    for line in lines:
        # Split the line by whitespace and take the first 3 elements
        parts = line.split()
        if len(parts) >= 3:  # Ensure the line has at least 3 parts
            # Convert the first three parts to a tuple and append to the list
            coords = (parts[0], parts[1], parts[2])
            coordinates.append(coords)
    return coordinates



def find_and_read_lines_from_string(content, target_line, num_of_atoms):
    """
    Finds a target line within the content string, and then reads n lines after the target line.

    :param content: String containing the content of the file.
    :param target_line: The line to search for in the content.
    :param n_lines_after: The number of lines to read after finding the target line.
    :return: A list of lines including the target line and the next n lines after it.
    """
    # Initialize an empty list to store the lines found after the target line
    lines_found = []

    # Split the content string into a list of lines for processing
    lines = content.split('\n')

    # Iterate through each line to find the target line
    for i, line in enumerate(lines):
        if target_line in line:
            # Once the target line is found, calculate the range of lines to extract
            start_index = i
            end_index = min(i + num_of_atoms + 1, len(lines))
            # Extract the specified range of lines, including the target line
            lines_found = lines[start_index+1:end_index]
            break  # Exit the loop since the target line has been found

    # Check if any lines were found after the target line
    if not lines_found:
        print(f"Target line '{target_line}' not found in the content.")
        return []


    return lines_found

def extract_coordinates_by_element(lines):
    coordinates_by_element = {}
    for line in lines:
        parts = line.split()
        if len(parts) >= 6:  # Ensure the line has enough parts to include coordinates and element symbol
            element_symbol = parts[-1]  # The element symbol is the last part of the line
            coords = (parts[0], parts[1], parts[2], parts[3], parts[4], parts[5])  # The coordinates are the first three parts
            # Add the coordinates to the list for this element in the dictionary
            if element_symbol in coordinates_by_element:
                coordinates_by_element[element_symbol].append(coords)
            else:
                coordinates_by_element[element_symbol] = [coords]
    return coordinates_by_element


    

if __name__ == '__main__':
    main()


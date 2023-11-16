__all__ = ['find_string']



def find_string(file, string):
    with open(file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if string in line:
                return i + 1
    return None



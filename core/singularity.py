__all__ = ['Singularity']

import os, json
import Arch 


class Singularity(object):
    def __init__(self, path_json):
        self.path_json = path_json 
        try:
            with open(self.path_json, 'r') as file_input:
                self.param = json.load(file_input)
        except:
            print("No such file or directory: '%s'"%self.path_json)
        else:
            print("The JSON file is imported successfully!")
            
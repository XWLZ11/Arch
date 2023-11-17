__all__ = ['Singularity']


import json
from ..lib.check import check_func

class Singularity(object):
    def __init__(self, path_json):
        self.path_json = path_json
        self.results = {} 
        try:
            with open(self.path_json, 'r') as file_input:
                self.param = json.load(file_input)
        except:
            print("No such file or directory: '%s'"%self.path_json)
        else:
            print("The JSON file is imported successfully!")
    
    def run(self):
        list_module = _find_inner_dict(self.param)
        for module in list_module:
            parts = module.split('.')
            name_module = '.'.join(parts[:-2])
            name_func = parts[-2]
            func = check_func(name_module, name_func)
            self.results.update({name_func: func(**self.param)})
        return self.results

def _find_inner_dict(dictionary, path="", result=None):
    if result is None:
        result = []
    if isinstance(dictionary, dict):
        for key, value in dictionary.items():
            if isinstance(value, dict):
                result = _find_inner_dict(value, path + key + '.', result)
            elif key == "switch" and value is True:
                result.append(path)
    return result


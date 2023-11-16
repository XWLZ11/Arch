__all__ = ['check_module', 'check_func']


import importlib, sys


def check_module(name_module):
    '''
    通用函数，用于校验是否将第三方库成功导入并加入命名空间。

    name_module:第三方库的完整名，可以包括其下的某个模块。

    举例：
    os = check_module('os')
    '''
    if name_module not in sys.modules:
        try:
            alias_module = importlib.import_module(name_module)
        except ModuleNotFoundError as e:
            print("Module name error, unable to import!")
            sys.exit(e)

    else:
        alias_module = sys.modules[name_module] 

    return alias_module


def check_func(name_module, name_func):
    '''
    通用函数，用于校验是否将第三方库的某一函数/方法成功导入并加入命名空间。

    name_module:第三方库的完整名，可以包括其下的某个模块。
    name_func:目标函数/方法的完整名

    举例：
    read = check_func('ase.io', 'read')
    '''    
    alias_module = check_module(name_module)
    try:
        if hasattr(alias_module, name_func):
            alias_func = getattr(alias_module, name_func)
        else:
            raise SystemExit(1)
    except SystemExit as e:
        print("Module or function name error, unable to import!")
        sys.exit(e)
        
    return alias_func


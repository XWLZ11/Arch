__all__ = ['unzipfile', 'zipfile']



from ...check import check_module, check_func

def unzipfile(dir_init="", dir_target=""):
    '''
    通用函数，用于解压zip格式文件。

    dir_init：初始目录路径
    dir_target：目标目录路径
    '''
    zf = check_module('zipfile')

    file_zip = zf.ZipFile(dir_init)
    file_unzip = file_zip.extractall(dir_target)
    file_zip.close()
    print("unzipfile all done!")

def zipfile(dir_init="", dir_target=""):
    '''
    通用函数，用于解压zip格式文件。

    dir_init：初始目录路径
    dir_target：目标目录路径
    '''
    zf = check_module('zipfile')
    os = check_module('os')
    
    filename_zip = os.path.dirname(dir_init) + '.zip'
    file_zip = zf.ZipFile(filename_zip, 'w', zf.ZIP_DEFLATED)

    for dirname, subfolders, filenames in os.walk(dir_init):
        for filename in filenames:
            # 构建文件的完整路径
            path_file = os.path.join(dirname, filename)
            # 构建解压后的路径
            path_extract = os.path.join(dir_target, os.path.relpath(path_file, dir_init))  
            # 创建文件夹结构
            os.makedirs(os.path.dirname(path_extract), exist_ok=True)
            file_zip.write(path_file, arcname=os.path.relpath(path_file, dir_init))

    file_zip.close()
    print("zipfile all done!")

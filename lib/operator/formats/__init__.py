'''
作用：不同坐标文件格式转换，文件压缩/解压

导入方式：
from Arch.lib.operator.formats import formats
'''

__all__ = ['pdb2xyz', 'pdb2lammpsdata', 'lammpstrj2lammpsdata', 'npy2poscar', 
           'poscar2xyz', 'xyz2npy', 'xyz2lammpstrj', 'lammpstrj2xyz', 'xsd2pdb', 
           'unzipfile', 'zipfile']


from ....core.singularity import Singularity
from .transform import *
from .zip import *

class formats(Singularity):
    '''
    通用类，用于处理不同格式的文件。
    '''
    def __init__(self, path_json):
        '''
        内部初始化方法，用于设置缺省值。
        '''
        super().__init__(path_json)
        self.dir_init = self.param.get("dir_init", "")
        self.dir_target = self.param.get("dir_target", "")

    def _checkerror(self, dir_init, dir_target):
        '''
        内部调用方法，用于核对目录路径是否正确。

        dir_init：初始目录路径
        dir_target：目标目录路径
        '''
        try:
            if dir_init == "" and self.dir_init != "":
                dir_init = self.dir_init
                dir_target = self.dir_target
            elif dir_init == "" and self.dir_init == "":
                raise SystemExit(1)
        except SystemExit as e:
            print("File or directory path error!")
            import sys
            sys.exit(e)
        else:
            if dir_target == "" and self.dir_target == "":
                dir_target = dir_init
            elif dir_target == "" and self.dir_target != "":
                dir_target = self.dir_target
            return dir_init, dir_target

    def pdb2xyz(self, dir_init="", dir_target="", names=['O', 'H']):
        '''
        通用方法，用于将pdb格式文件转换为xyz格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径
        names - 原子名列表（赝势顺序）

        用法举例：
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        if names == ['O', 'H']:
            names = self.param.get("names_pdb2xyz", ['O', 'H'])
        pdb2xyz(dir_init, dir_target, names)

    def pdb2lammpsdata(self, dir_init="", dir_target="", names=['O', 'H']):
        '''
        通用方法，用于将pdb格式文件转换为lammpsdata格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径
        names - 原子名列表（赝势顺序）

        用法举例：
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        if names == ['O', 'H']:
            names = self.param.get("names_pdb2lammpsdata", ['O', 'H'])
        pdb2lammpsdata(dir_init, dir_target, names)

    def lammpstrj2lammpsdata(self, dir_init="", dir_target=""):
        '''
        通用方法，用于将pdb格式文件转换为lammpsdata格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径

        用法举例：
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        lammpstrj2lammpsdata(dir_init, dir_target)


    def npy2poscar(self, dir_init="", dir_target="", multi=False, record=0):
        '''
        通用方法，用于将npy格式文件转换为POSCAR格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径
        multi：是否批量转换npy格式
        record：用于记录批量转换npy格式的进度

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
    
        list_dir = os.listdir(S.param["dir_init"])
        file_pos = formats(path_json)
        record = 0
        for dirname in list_dir:
            record = file_pos.npy2poscar(S.param["dir_init"]+dirname, multi=True, record=record)
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        if multi is False:
            multi = self.param.get("multi_npy2poscar", False)
        if record == 0:
            record = self.param.get("record_npy2poscar", 0)
        record = npy2poscar(dir_init, dir_target, multi, record)
        return record

    def poscar2xyz(self, dir_init="", dir_target=""):
        '''
        通用方法，用于将POSCAR格式文件转换为xyz格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.poscar2xyz()
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        poscar2xyz(dir_init, dir_target)

    def xyz2npy(self, dir_init="", dir_target="", multi=False, record=0):
        '''
        通用方法，用于将xyz格式文件转换为npy格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.xyz2npy()
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        if multi is False:
            multi = self.param.get("multi_xyz2npy", False)
        if record == 0:
            record = self.param.get("record_xyz2npy", 0)
        xyz2npy(dir_init, dir_target, multi, record)
    
    def xyz2lammpstrj(self, dir_init="", dir_target=""):
        '''
        通用方法，用于将xyz格式文件转换为lammpstrj格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.xyz2lammpstrj()
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        xyz2lammpstrj(dir_init, dir_target)

    def lammpstrj2xyz(self, dir_init="", dir_target="", names=['O', 'H']):
        '''
        通用方法，用于将xyz格式文件转换为lammpstrj格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.lammpstrj2xyz()
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        if names == ['O', 'H']:
            names = self.param.get("names_lammpstrj2xyz", ['O', 'H'])
        lammpstrj2xyz(dir_init, dir_target, names)

    def xsd2pdb(self, dir_init="", dir_target="", filename="system.xsd"):
        '''
        通用方法，用于将xsd格式文件转换为pdb格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.xsd2pdb()
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)        
        if filename == 'system.xsd':
            filename = self.param.get("filename_xsd2pdb", "system.xsd")
        xsd2pdb(dir_init, dir_target, filename)    

    def unzipfile(self, dir_init="", dir_target=""):
        '''
        通用函数，用于解压zip格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.unzipfile()        
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        unzipfile(dir_init, dir_target)

    def zipfile(self, dir_init="", dir_target=""):
        '''
        通用函数，用于解压zip格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.zipfile()
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        zipfile(dir_init, dir_target)
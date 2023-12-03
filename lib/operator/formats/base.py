from ....core.base import Singularity
from .transform import *
from .zip import *

class Formats(Singularity):
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

    def pdb2xyz(self, **kwargs):
        '''
        通用方法，用于将pdb格式文件转换为xyz格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径
        names - 原子名列表（赝势顺序）

        用法举例：
        '''
        pdb2xyz(**self.param)

    def pdb2lammpsdata(self, **kwargs):
        '''
        通用方法，用于将pdb格式文件转换为lammpsdata格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径
        names - 原子名列表（赝势顺序）

        用法举例：
        '''
        pdb2lammpsdata(**self.param)

    def lammpstrj2lammpsdata(self, **kwargs):
        '''
        通用方法，用于将pdb格式文件转换为lammpsdata格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径

        用法举例：
        '''
        lammpstrj2lammpsdata(**self.param)


    def npy2poscar(self, **kwargs):
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
        record = npy2poscar(**self.param)
        return record

    def poscar2xyz(self, **kwargs):
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
        poscar2xyz(**self.param)

    def qelog2npy(self, **kwargs):
        '''
        通用方法，用于将xyz格式文件转换为npy格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.qelog2npy()
        '''
        qelog2npy(**self.param)

    def xyz2npy(self, **kwargs):
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
        xyz2npy(**self.param)
    
    def xyz2lammpstrj(self, **kwargs):
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
        xyz2lammpstrj(**self.param)

    def lammpstrj2xyz(self, **kwargs):
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
        lammpstrj2xyz(**self.param)

    def xsd2pdb(self, **kwargs):
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

        xsd2pdb(**self.param)    
    def lammpstrj2lammpstrj(self, **kwargs):
        '''
        通用方法，用于将xsd格式文件转换为pdb格式文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        file_pos = formats(path_json)
        file_pos.lammpstrj2lammpstrj()
        '''
        lammpstrj2lammpstrj(**self.param)        

    def unzipfile(self, **kwargs):
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
        #dir_init, dir_target = self._checkerror(dir_init, dir_target)
        unzipfile(**self.param)

    def zipfile(self, **kwargs):
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
        #dir_init, dir_target = self._checkerror(dir_init, dir_target)
        zipfile(**self.param)
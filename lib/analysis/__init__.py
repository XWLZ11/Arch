

__all__ = ['compute_RDF', 'find_string']


from ...core.singularity import Singularity
from .MD import *
from .seek import *

class MD(Singularity):
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

    def compute_RDF(self, dir_init="", dir_target="", file_lammpstrj='dump.lammpstrj', 
                    selection1='type 1', selection2='type 2', start=0, stop=-1, step=1):
        '''
        通用方法，用于将pdb格式文件转换为xyz格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径


        用法举例：
        '''
        dir_init, dir_target = self._checkerror(dir_init, dir_target)
        if file_lammpstrj == 'dump.lammpstrj':
            file_lammpstrj = self.param.get("name_file_lammpstrj", "dump.lammpstrj")
        if selection1 == 'type 1':
            selection1 = self.param.get("selection1", "type 1")
        if selection2 == 'type 2':
            selection2 = self.param.get("selection2", "type 2")
        if start == 0:
            start = self.param.get("start_compute_RDF", 0)
        if stop == -1:
            start = self.param.get("stop_compute_RDF", -1)
        if step == 1:
            start = self.param.get("step_compute_RDF", 0)        
        [r, g_r] = compute_RDF(dir_init, dir_target, file_lammpstrj, selection1, 
                                selection2, start, stop, step)
        
        return [r, g_r]

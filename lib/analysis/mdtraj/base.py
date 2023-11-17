
from ....core.base import Singularity
from .RDF import *

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

    def rdf(self, **kwargs):
        '''
        通用方法，用于将pdb格式文件转换为xyz格式文件。

        dir_init - 初始目录路径
        dir_target - 目标目录路径


        用法举例：
        '''
        [r, g_r] = rdf(**self.param)

        return [r, g_r]

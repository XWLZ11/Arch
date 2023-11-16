'''
导入方式：
from Arch.lib.operator.splice import splice
'''

__all__ = ['splice_input']

from ....core.singularity import Singularity
from ...check import check_module, check_func
from ...analysis.seek import find_string


class splice(Singularity):
    '''
    通用类，用于拼接不同文件。
    '''
    def __init__(self, path_json):
        '''
        内部初始化方法，用于设置缺省值。
        '''
        super().__init__(path_json)
        self.dir_init = self.param.get("dir_init", "")
        self.dir_target = self.param.get("dir_target", "")
        self.file_input = self.param.get("file_input", "")
        
    def _checkerror(self, dir_init, dir_target, file_input):
        '''
        内部调用方法，用于核对目录路径是否正确。

        dir_init：初始目录路径
        dir_target：目标目录路径
        '''
        try:
            if dir_init == "" and self.dir_init != "" or file_input == "" and self.file_input != "":
                dir_init = self.dir_init
                dir_target = self.dir_target
                file_input = self.file_input
            elif dir_init == "" and self.dir_init == "" or file_input == "" and self.file_input == "":
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
            return dir_init, dir_target, file_input

        
    def splice_input(self, dir_init="", dir_target="", file_input=""):
        '''
        通用方法，用于将两个文件拼接为一个输入文件。

        dir_init：初始目录路径
        dir_target：目标目录路径

        用法举例：
        path_json = '/home/tsingularity/Group/scripts/my_scripts/para_input.json'
        S = Arch.Singularity(path_json)
        T = splice(path_json)
        T.splice_input()
        '''
        dir_init, dir_target, file_input = self._checkerror(dir_init, dir_target, file_input)
        os = check_module('os')

        try:
            os.makedirs(dir_target+'input/')
        except FileExistsError as e:
            import sys 
            sys.exit(e)
        else:
            list_filenames = os.listdir(dir_init+'QE_xyz/')
            line_number_input = find_string(file_input, "ATOMIC_POSITIONS")

            for filename in list_filenames:
                with open(file_input, 'r') as input_file, open(dir_target+'input/'+filename+'.in', 'w') as output_file:
                    lines = input_file.readlines()[:line_number_input]
                    output_file.writelines(lines)
                with open(dir_init+'QE_xyz/'+filename, 'r') as input_file, open(dir_target+'input/'+filename+'.in', 'a') as output_file:
                    lines = input_file.readlines()[2:]
                    output_file.writelines(lines)
            print("splice all done!")
# The content of build f2py modules .py
  
from numpy.distutils.core import Extension
ext = Extension (name = "loss_gain_f2py",
                 sources = ["Loss_Gain.f90"], 
                 extra_compile_args=["-O3",
                                     "-ffast-math",
                                     "-fopenmp"],
                 extra_link_args=["-lgomp",
                                  "-static",
                                  "-llapack",
                                  "-lblas"])
if __name__ == '__main__':
    from numpy.distutils.core import setup 
    setup (name = "loss_gain_f2py" ,
           ext_modules = [ext])

    
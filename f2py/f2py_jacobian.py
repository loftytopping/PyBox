# The content of build f2py modules .py
  
from numpy.distutils.core import Extension
ext = Extension (name = "jacobian_f2py",
                 sources = ["Jacobian.f90"], 
                 extra_compile_args=["-O3",
                                     "-ffast-math",
                                     "-fopenmp"],
                 extra_link_args=["-lgomp",
                                  "-static",
                                  "-llapack",
                                  "-lblas"])
if __name__ == '__main__':
    from numpy.distutils.core import setup 
    setup (name = "jacobian_f2py" ,
           ext_modules = [ext])

    
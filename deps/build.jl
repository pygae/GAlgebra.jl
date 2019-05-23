# using PyCall;

# if !haskey(ENV, "USE_CUSTEM_GALGEBRA")
#   run(PyCall.python_cmd(`-m pip install sympy==1.3 -e git+https://github.com/pygae/galgebra.git#egg=galgebra`))
# end

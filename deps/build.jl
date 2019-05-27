using PyCall;

try
  pyimport("galgebra")
catch e
  run(PyCall.python_cmd(`-m pip install -e git+https://github.com/pygae/galgebra.git --user`))
  run(PyCall.python_cmd(`-m pip install sympy==1.3 --user`))
end
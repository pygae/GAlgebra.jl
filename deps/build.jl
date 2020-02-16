using PyCall;

try
  pyimport("galgebra")
catch e
  run(PyCall.python_cmd(`-m pip install galgebra --user`))
end
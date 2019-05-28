macro define_show(type)
  @eval begin
      Base.show(io::IO, x::$type) = print(io, pystr(x.o))
      Base.show(io::IO, ::MIME"text/plain", x::$type) = print(io, pystr(x.o))
      Base.show(io::IO, ::MIME"text/latex", x::$type) = print(io, "\\begin{align*}" * galgebra.printer.latex(x.o) * "\\end{align*}")
  end
end

macro delegate_properties(type, obj_field)
  @eval begin
      Base.convert(::Type{$type}, o::PyCall.PyObject) = $type(o)
      PyCall.PyObject(o::$type) = PyCall.PyObject(getfield(o, $obj_field))

      function Base.getproperty(o::$type, s::AbstractString)
          if s == String($obj_field)
              return getfield(o, $obj_field)
          else
              return getproperty(getfield(o, $obj_field), s)
          end
      end
      
      Base.getproperty(o::$type, s::Symbol) = getproperty(o, String(s))
      
      Base.propertynames(o::$type) = map(x->Symbol(first(x)),
                                      pycall(inspect."getmembers", PyObject, getfield(o, $obj_field)))
      
      # avoiding method ambiguity
      Base.setproperty!(o::$type, s::Symbol, v) = _setproperty!(o,s,v)
      Base.setproperty!(o::$type, s::AbstractString, v) = _setproperty!(o,s,v)
      
      function _setproperty!(o::$type, s::Union{Symbol,AbstractString}, v)
          obj = getfield(o, $obj_field)
          setproperty!(obj, s, v)
          o
      end

      hasproperty(o::$type, s::Symbol) = hasproperty(getfield(o, $obj_field), s)
      hasproperty(o::$type, s::AbstractString) = hasproperty(getfield(o, $obj_field), s)
  end
end

macro delegate_doc(type)
  @eval begin
      # Expose Python docstrings to the Julia doc system
      Docs.getdoc(x::$type) = Text(convert(String, x."__doc__"))
      Docs.Binding(x::$type, s::Symbol) = getproperty(x, s)
  end
end

macro define_op(type, op, method)
  @eval begin
      @pure $op(x::$type, y::$type) = x.$method(y)
  end
end

macro define_lop(type, rtype, op, lmethod)
  @eval begin
      @pure $op(x::$type, y::$rtype) = x.$lmethod(y)
  end
end
              
macro define_rop(type, ltype, op, rmethod)
  @eval begin
      @pure $op(x::$ltype, y::$type) = y.$rmethod(x)
  end
end

macro define_unary_op(type, op, method)
  @eval begin
      @pure $op(x::$type) = x.$method()
  end
end

macro define_postfix_symbol(super_script)
  @eval begin
      const $super_script = () -> String(:super_script)
      export $super_script
  end
end

macro define_postfix_op(type, super_script, func)
  @eval begin
      Base.:(*)(x::$type,::typeof($super_script)) = $func(x)
      Base.:(^)(x::$type,::typeof($super_script)) = $func(x)
  end
end

macro define_show(type)
  return quote
      Base.show(io::IO, x::$(esc(type))) = print(io, pystr(x.o))
      Base.show(io::IO, ::MIME"text/plain", x::$(esc(type))) = print(io, pystr(x.o))
      Base.show(io::IO, ::MIME"text/latex", x::$(esc(type))) = print(io, "\\begin{align*}" * galgebra.printer.latex(x.o) * "\\end{align*}")
  end
end

macro delegate_properties(type, obj_field)
  return quote
      Base.convert(::Type{$(esc(type))}, o::PyCall.PyObject) = $(esc(type))(o)
      PyCall.PyObject(o::$(esc(type))) = PyCall.PyObject(getfield(o, $obj_field))

      function Base.getproperty(o::$(esc(type)), s::AbstractString)
          if s == String($obj_field)
              return getfield(o, $obj_field)
          else
              return getproperty(getfield(o, $obj_field), s)
          end
      end
      
      Base.getproperty(o::$(esc(type)), s::Symbol) = getproperty(o, String(s))
      
      Base.propertynames(o::$(esc(type))) = map(x->Symbol(first(x)),
                                      pycall(inspect."getmembers", PyObject, getfield(o, $obj_field)))
      
      # avoiding method ambiguity
      Base.setproperty!(o::$(esc(type)), s::Symbol, v) = _setproperty!(o,s,v)
      Base.setproperty!(o::$(esc(type)), s::AbstractString, v) = _setproperty!(o,s,v)
      
      function _setproperty!(o::$(esc(type)), s::Union{Symbol,AbstractString}, v)
          obj = getfield(o, $obj_field)
          setproperty!(obj, s, v)
          o
      end

      hasproperty(o::$(esc(type)), s::Symbol) = hasproperty(getfield(o, $obj_field), s)
      hasproperty(o::$(esc(type)), s::AbstractString) = hasproperty(getfield(o, $obj_field), s)
  end
end

macro delegate_doc(type)
  return quote
      # Expose Python docstrings to the Julia doc system
      Docs.getdoc(x::$(esc(type))) = Text(convert(String, x."__doc__"))
      Docs.Binding(x::$(esc(type)), s::Symbol) = getproperty(x, s)
  end
end

macro define_op(type, op, method)
  return quote
    Core.@__doc__ @pure $(esc(op))(x::$(esc(type)), y::$(esc(type))) = x.$method(y)
  end
end

macro define_lop(type, rtype, op, lmethod)
  return quote
    Core.@__doc__ @pure $(esc(op))(x::$(esc(type)), y::$(esc(rtype))) = x.$lmethod(y)
  end
end
              
macro define_rop(type, ltype, op, rmethod)
  return quote
    Core.@__doc__ @pure $(esc(op))(x::$(esc(ltype)), y::$(esc(type))) = y.$rmethod(x)
  end
end

macro define_unary_op(type, op, method)
  return quote
    Core.@__doc__ @pure $(esc(op))(x::$(esc(type))) = x.$method()
  end
end

macro define_postfix_symbol(super_script)
  return quote
      const $(esc(super_script)) = () -> String(:super_script)
      export $(esc(super_script))
  end
end

macro define_postfix_op(type, super_script, func)
  return quote
      Base.:(*)(x::$(esc(type)),::typeof($(esc(super_script)))) = $func(x)
      Base.:(^)(x::$(esc(type)),::typeof($(esc(super_script)))) = $func(x)
  end
end

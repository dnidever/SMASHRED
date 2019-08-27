function make_structure,str1,num
  
schema = str1[0]
struct_assign,{dum:''},schema
return,replicate(schema,num)

end

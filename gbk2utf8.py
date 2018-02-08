import sys
file_name=sys.argv[1]
f=open(file_name,'rb')
st=f.read().decode('gbk')
f=open(file_name,'wb')
f.write(st.encode('utf-8'))
f.close()


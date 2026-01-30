content = "fastapi\nuvicorn\npydantic\nsqlalchemy\n"
with open('cf-requirements.txt','w',encoding='utf-8') as f:
    f.write(content)
print('Wrote cf-requirements.txt')

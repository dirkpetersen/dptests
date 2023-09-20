from litestar import Litestar, get

# look at these 
# https://www.infoworld.com/article/3700689/3-python-web-frameworks-for-beautiful-front-ends.html
# https://github.com/litestar-org/litestar-fullstack
# https://github.com/v3ss0n/awesome-starlite/

# https://anvil.works/
# or Reflex.dev, NiceGUI 

@get("/")
async def hello_world() -> str:
    return "Hello, world!"

app = Litestar([hello_world])

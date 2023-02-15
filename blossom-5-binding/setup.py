from cmaketools import setup

setup(
    name="blossom_5_binding",
    version="0.1.0",
    author="Alexis Shaw",
    author_email="alexis@alexisshaw.com",
    description="A python binding for Blossom-V using pybind11",
    url="http://alexisshaw.com",
    license="Proprietary / MIT",
    src_dir="src",
    ext_module_hing=r"pybind11_add_module",
    has_package_data=False,
)

# setup(
#     name="pybind_example",
#     version="0.1.0",
#     author="Takeshi (Kesh) Ikuma",
#     author_email="tikuma@gmail.com",
#     description="A test package using pybind11",
#     url="https://github.com/python-cmaketools/pybind-example.git",
#     license="MIT License",
#     src_dir="src",
#     ext_module_hint=r"pybind11_add_module",
#     has_package_data=False,
# )

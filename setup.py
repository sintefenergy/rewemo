from distutils.core import setup
setup(
	name = "rewemo",
	version = "0.1.0",
	description = "Renewable energy time series from numerical weather model data",
	author = "Harald G Svendsen",
	author_email = "harald.svendsen@sintef.no",
	license = "MIT License (http://opensource.org/licenses/MIT)",
	classifiers = [
		"Programming Language :: Python :: 3",
		"Development Status :: 4 - Beta",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License"
		],
#	packages = ["rewemo"],
	package_dir= {"": "src/rewemo"},
	python_requires = ">=3.7"
)
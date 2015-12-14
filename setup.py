from setuptools import setup

setup(
    name='rnacount',
    description='Count multiplexed shRNA reads',
    author='Adrian Seyboldt',
    author_email='adrian.seyboldt@gmail.com',
    url='https://github.com/qbicsoftware/rnacount',
    scripts = ['rnacount.py'],
    install_requires = [
        'numpy',
        'pandas',
    ],
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
)

from setuptools import setup, find_packages

setup(
    name="q2-periodontitis",
    version="2026.1.0",
    packages=find_packages(),
    author="Fabiano Menegidio",
    author_email="fabiano.menegidio@example.com",
    description="QIIME 2 plugin for subgingival biofilm ecological characterization, Socransky complexes collapse, and dysbiosis indices.",
    license="MIT",
    url="https://github.com/LaBiOmicS/q2-periodontitis",
    entry_points={
        "qiime2.plugins": [
            "q2-periodontitis=q2_periodontitis.plugin_setup:plugin"
        ]
    },
    package_data={
        "q2_periodontitis": [
            "citations.bib",
            "data/*",
            "assets/web/*",
            "assets/web/css/*",
            "assets/web/js/*"
        ]
    },
    zip_safe=False,
)

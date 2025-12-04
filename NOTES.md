Build:
```
python -m pip install --upgrade build
python -m build
```

Install local:
```
pip install dist/nfixplanet-0.1.0-py3-none-any.whl
```

Once test passes:
```
pip install --upgrade twine
twine upload dist/*
```

OR use uv to do everything
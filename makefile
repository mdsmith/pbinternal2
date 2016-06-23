auto-pep8:
	find bin -name "*.py" -exec autopep8 -i --ignore=E501,E265,E731,E402,W292 {} \;

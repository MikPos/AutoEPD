import importlib

tests = ["test_book",
	# "test_daniel",
	"test_formose", "test_long_klaus"]

exs = []
for t in tests:
	m = importlib.import_module(t)
	exs.extend(m.exAll)

for E in exs:
	E().run()

import importlib

def import_eos(eos):
    importlib.invalidate_caches()
    try:
        importlib.import_module('ComPASS.eos.%s' % eos)
    except Exception as e:
        print(e)
        return False
    return True

def test_import_all_eos():
    assert import_eos('water2ph')
    assert import_eos('diphasic')

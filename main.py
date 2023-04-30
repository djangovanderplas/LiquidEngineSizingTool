from InitialSizing import engine_sizer


test_engine = engine_sizer.Engine('config.yaml')
print(test_engine.isp)
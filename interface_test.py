import ALEAF_interface as ALI

ALEAF_path = "/home/biegelk/kb_aleaf/data/LC_GEP/ERCOT/ALEAF_ERCOT.xlsx"
new_ALEAF_path = "/home/biegelk/kb_aleaf/data/LC_GEP/ERCOT/ALEAF_ERCOT_updated.xlsx"

sys_port = ALI.load_ALEAF_system_portfolio(ALEAF_path)
non_port = ALI.load_non_portfolio_sheets(ALEAF_path)

print(sys_port)
print("\n")
print(non_port)

ALI.overwrite_ALEAF_system_file(new_ALEAF_path, sys_port, non_port)

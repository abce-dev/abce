import ABCEfunctions as af
import yaml

settings = yaml.load(open("/home/biegelk/abce/settings.yml", "r"), Loader=yaml.FullLoader)
output_dir = "/home/biegelk/abce/outputs/ABCE_base"

af.process_outputs(settings, output_dir)


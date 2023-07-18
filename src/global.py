import logging


# Configure the root logger
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Create a StreamHandler and set its output to sys.stdout
console_handler = logging.StreamHandler(sys.stdout)

# Optionally, you can set the logging level for the handler
console_handler.setLevel(logging.DEBUG)

# Add the handler to the root logger
logging.getLogger().addHandler(console_handler)
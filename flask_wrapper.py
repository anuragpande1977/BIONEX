from flask import Flask, redirect

app = Flask(__name__)

@app.route("/")
def index():
    # Redirect to the Streamlit app URL
    return redirect("https://bionex.streamlit.app")  # Replace with your Streamlit app URL

import os

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))  # Get the port from the environment or use 5000 as default
    app.run(debug=True, host="0.0.0.0", port=port)


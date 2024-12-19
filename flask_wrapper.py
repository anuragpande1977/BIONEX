from flask import Flask, redirect

app = Flask(__name__)

@app.route("/")
def index():
    # Redirect directly to Streamlit app
    return redirect("https://bionex.streamlit.app")

import os

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(debug=True, host="0.0.0.0", port=port)


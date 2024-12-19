from flask import Flask, render_template_string

app = Flask(__name__)

@app.route("/")
def index():
    # Google Tag script for <head>
    gtag_script = """
    <!-- Google tag (gtag.js) -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-93XG1G92HZ"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-93XG1G92HZ');
    </script>
    """

    # Simple HTML template
    return render_template_string(f"""
    <!DOCTYPE html>
    <html>
    <head>
        {gtag_script}
        <title>Flask App with Google Tag</title>
    </head>
    <body>
        <h1>Google Tag Test</h1>
        <p>If you can see this, the Google Tag has been injected into the page.</p>
    </body>
    </html>
    """)

import os

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))  # Use PORT from Render or default to 5000
    app.run(debug=True, host="0.0.0.0", port=port)

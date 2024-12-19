from flask import Flask, render_template_string

app = Flask(__name__)

@app.route("/")
def index():
    # Google Tag Manager (GTM) script for <head>
    gtm_head = """
    <!-- Google Tag Manager -->
    <script>
      (function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':
      new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],
      j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src=
      'https://www.googletagmanager.com/gtm.js?id=GTM-PZKDTSN8';f.parentNode.insertBefore(j,f);
      })(window,document,'script','dataLayer','GTM-PZKDTSN8');
    </script>
    <!-- End Google Tag Manager -->
    """

    # Google Tag Manager (GTM) <noscript> tag for <body>
    gtm_body = """
    <!-- Google Tag Manager (noscript) -->
    <noscript>
      <iframe src="https://www.googletagmanager.com/ns.html?id=GTM-PZKDTSN8"
      height="0" width="0" style="display:none;visibility:hidden"></iframe>
    </noscript>
    <!-- End Google Tag Manager (noscript) -->
    """

    # Streamlit app URL (update with your Streamlit app's public URL)
    streamlit_app_url = "https://bionex.streamlit.app"  # Replace with your Streamlit app URL

    # HTML template with GTM and embedded Streamlit app
    return render_template_string(f"""
    <!DOCTYPE html>
    <html>
    <head>
        {gtm_head}
    </head>
    <body>
        {gtm_body}
        <iframe src="{streamlit_app_url}" width="100%" height="100%" style="border:none;"></iframe>
    </body>
    </html>
    """)

if __name__ == "__main__":
    app.run(debug=True, port=5000)

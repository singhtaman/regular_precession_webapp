"""Entry point for the redesigned Regular Precession explorer.

Vercel loads the WSGI callable named ``app`` from this file. Locally, run
``python app.py`` and open http://127.0.0.1:8050.

The Render deployment is separate and still serves RP_map_all_app.py from `main`.
"""
import os

from webapp.main import app as dash_app, server

app = server  # WSGI application (Flask) for Vercel and other WSGI servers

if __name__ == "__main__":
    dash_app.run(host=os.environ.get("HOST", "127.0.0.1"),
                 port=int(os.environ.get("PORT", 8050)),
                 debug=os.environ.get("DASH_DEBUG") == "1")

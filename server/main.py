from flask import Flask, render_template

app = Flask(__name__)

@app.route("/")
def base_get():
        return render_template("notebook.html")

if __name__ == '__main__':
        app.run(debug=False, host="0.0.0.0", port=5000)

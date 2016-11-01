from datetime import timedelta

from json import dumps, loads
from functools import update_wrapper
from flask import Flask, jsonify, request, make_response, send_file, current_app
from werkzeug.contrib.fixers import ProxyFix

from main import process, process_batch_folding
from plotter import RNAPlotter


app = Flask(__name__)
app.wsgi_app = ProxyFix(app.wsgi_app)
app.config['PROPAGATE_EXCEPTIONS'] = True


def crossdomain(origin=None, methods=None, headers=None, 
	max_age=21600, attach_to_all=True, automatic_options=True):
    if methods is not None:
        methods = ', '.join(sorted(x.upper() for x in methods))
    if headers is not None and not isinstance(headers, basestring):
        headers = ', '.join(x.upper() for x in headers)
    if not isinstance(origin, basestring):
        origin = ', '.join(origin)
    if isinstance(max_age, timedelta):
        max_age = max_age.total_seconds()

    def get_methods():
        if methods is not None:
            return methods
        options_resp = current_app.make_default_options_response()
        return options_resp.headers['allow']

    def decorator(f):
        def wrapped_function(*args, **kwargs):
            if automatic_options and request.method == 'OPTIONS':
                resp = current_app.make_default_options_response()
            else:
                resp = make_response(f(*args, **kwargs))
            if not attach_to_all and request.method != 'OPTIONS':
                return resp

            h = resp.headers

            h['Access-Control-Allow-Origin'] = origin
            h['Access-Control-Allow-Methods'] = get_methods()
            h['Access-Control-Max-Age'] = str(max_age)
            if headers is not None:
                h['Access-Control-Allow-Headers'] = headers
            return resp

        f.provide_automatic_options = False
        return update_wrapper(wrapped_function, f)

    return decorator


@app.route('/api/process', methods=['POST'])
@crossdomain(origin='*')
def process_post_request():
	params = request.get_json()
	return jsonify(process(params))


@app.route('/api/folding', methods=['POST'])
@crossdomain(origin='*')
def process_folding():
    params = request.get_json()
    return jsonify(process_batch_folding(params))


@app.route('/api/get_ss_image', methods=['POST'])
@crossdomain(origin='*')
def process_ss():
    params = request.get_json()
    plot_file_name = RNAPlotter.from_json(params)
    return send_file(plot_file_name, as_attachment=True, attachment_filename=plot_file_name)


if __name__ == '__main__':
	app.run(host='0.0.0.0', port=9000)

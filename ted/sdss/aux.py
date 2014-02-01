#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sun 29 Sep 2013
#   Initial build.
#

import os

def HTML_create_cutout_display(cutout_files, path=None):
    """Creates HTML document that displays the cutout images."""
    if path is None:
        return

    HTML_head = '''\
<!doctype html />
<html>
<head>
<title>Cutouts</title>
<style rel="stylesheet" type="text/css">
span {
    display: block;
    float: left;
    border: 1px #ddd solid;
    background: #f6f6f6;
    padding: 10px;
    margin: 10px 10px;
    text-align: center;
}
div.clear {
    clear: both;
}
span img {
    border: 5px #fff solid;
    margin: 0;
    padding: 0;
}
div {
    text-align: center;
}
div span img {
    width: auto;
}
pre {
    display:block;
}
</style>
</head>
<body>
'''
    HTML_foot = '''\
</body>
</html>
'''
    HTML_body = ''

    # Add cutout summary
    ifname_summary = os.path.join(path, 'summary.txt')
    if os.path.isfile(ifname_summary):
        print 'Found {}'.format(ifname_summary)
        with open(ifname_summary, 'r') as fsock:
            HTML_body += '\n<pre>{}</pre>'.format(fsock.read())

    # Add cutou statistics
    ifname_images = [
        os.path.join(path, 'stats.png'),
        # os.path.join(path, 'pixel_indices_joint.png'),
        os.path.join(path, 'coverage.png'),
        os.path.join(path, '../max_cutout_size.png'),
    ]
    for ifname_image in ifname_images:
        if os.path.isfile(ifname_image):
            print 'Found {}'.format(ifname_image)
            HTML_body += HTML_image(ifname_image, alt=None)

    # Add each cutout image
    for fname in cutout_files:
        HTML_body += '''

<span>
<img src="png/{0}" alt="{1}" title="{1}" /><br />
<em>{2}</em>
</span>'''.format(fname, os.path.splitext(fname)[0], fname.split('T')[0])

    # Combine HTML fragments and save the result
    HTML_content = HTML_head + HTML_body + HTML_foot
    ofname_html = os.path.join(path, 'index.html')
    with open(ofname_html, 'w+') as fsock:
        fsock.write(HTML_content)

    # Return True if file was created
    return os.path.isfile(ofname_html)


def HTML_image(ifname, alt=None):
    return '''

<div>
  <span><img src="{0}" alt="{1}" title="{1}" /></span>
</div>
<div class="clear"></div>'''.format(ifname, '' if alt is None else alt)


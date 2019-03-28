This API uses [Zappa](https://www.zappa.io/), a serverless python web framework built on top of Flask, for deployment to AWS.

## Deploying

```
$ virtualenv env
$ source env/bin/activate
$ pip install -r requirements.txt
$ zappa deploy
$ zappa update # further updates
```

## Testing

Local testing:

```
pytest
```

Using cURL to test Zappa deployed API endpoint:

```
$ curl -XGET -G -d "lat1=-37.57037203" -d "lon1=144.25295244" \ 
>               -d "azimuth1to2=306.520537" -d "ell_dist=54972.271" \
>               -d "from_angle_type=dms" -d "to_angle_type=dms" \
> https://<YOU_API_GATEWAY_ENDPOINT>/dev/vincdir
{"azimuth2to1":127.10250671432,"lat2":-37.3910156124268,"lon2":143.5535383883988}
```

import struct

def create_audio_event(audio_chunk):
    """Create an event stream message containing audio data."""
    headers = [
        (':content-type', 'application/octet-stream'),
        (':event-type', 'AudioEvent'),
        (':message-type', 'event')
    ]
    
    msg = build_event_message(headers, audio_chunk)
    return msg

def build_event_message(headers, payload):
    """Build event stream message with headers and payload."""
    # Calculate total header length
    header_section_length = sum(
        len(k) + len(str(v)) + 4  # +4 for header name/value length bytes
        for k, v in headers
    )
    
    # Calculate total message length
    total_length = (
        4 +  # prelude total bytes
        4 +  # prelude headers bytes
        4 +  # prelude CRC
        header_section_length +
        len(payload) +
        4  # message CRC
    )
    
    # Build message
    msg = bytearray()
    
    # Add prelude
    msg.extend(struct.pack('>I', total_length))
    msg.extend(struct.pack('>I', header_section_length))
    msg.extend(struct.pack('>I', calculate_crc32(msg[:8])))
    
    # Add headers
    for key, value in headers:
        msg.extend(struct.pack('B', len(key)))
        msg.extend(key.encode('utf-8'))
        msg.extend(struct.pack('B', 7))  # String type
        msg.extend(struct.pack('B', len(str(value))))
        msg.extend(value.encode('utf-8'))
    
    # Add payload
    msg.extend(payload)
    
    # Add message CRC
    msg.extend(struct.pack('>I', calculate_crc32(msg)))
    
    return bytes(msg)

def decode_event(msg):
    """Decode an event stream message into headers and payload."""
    pos = 0
    
    # Read prelude
    total_length = struct.unpack('>I', msg[pos:pos+4])[0]
    pos += 4
    headers_length = struct.unpack('>I', msg[pos:pos+4])[0]
    pos += 4
    prelude_crc = struct.unpack('>I', msg[pos:pos+4])[0]
    pos += 4
    
    # Read headers
    headers = {}
    headers_end = pos + headers_length
    while pos < headers_end:
        name_length = msg[pos]
        pos += 1
        name = msg[pos:pos+name_length].decode('utf-8')
        pos += name_length
        
        value_type = msg[pos]
        pos += 1
        value_length = msg[pos]
        pos += 1
        value = msg[pos:pos+value_length].decode('utf-8')
        pos += value_length
        
        headers[name] = value
    
    # Read payload
    payload = msg[pos:-4]  # Last 4 bytes are message CRC
    
    # If this is a JSON payload, decode it
    if headers.get(':content-type') == 'application/json':
        import json
        payload = json.loads(payload)
    
    return headers, payload

def calculate_crc32(data):
    """Calculate CRC32 checksum."""
    crc = 0xFFFFFFFF
    for b in data:
        crc ^= b
        for _ in range(8):
            crc = (crc >> 1) ^ 0xEDB88320 if crc & 1 else crc >> 1
    return crc ^ 0xFFFFFFFF

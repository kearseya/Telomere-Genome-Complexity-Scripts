import pytesseract
pytesseract.pytesseract.tesseract_cmd = r'tesseract'
print(pytesseract.image_to_string(r'large_test_crop.png'))
